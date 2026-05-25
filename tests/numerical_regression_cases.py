import json
import os
from pathlib import Path

for _name in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_name, "1")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/euler_elastica_matplotlib")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import numpy as np

from core import config
from core.analysis_solve import AnalysisFunc
from core.batch_cases import BATCH_TEST_CASES
from core.opensees_beam import get_fem_result
from core.RK import get_rk_result
from critical_boundary_fem import AdaptiveFemAnalyzer
from run_boundary_cpu_max import (
    _available_cpus,
    _build_tasks,
    _case_metadata,
    _choose_start_method,
    _group_rows,
    _parse_workers,
    _run_parallel_tasks,
)


NUMERIC_COLUMNS = ("s", "theta", "w", "x", "M", "V", "N")
BASELINE_ROOT = Path(__file__).resolve().parent / "numerical_baselines"

SINGLE_VAR_BOUNDARY_CASES = [1, 2, 4, 5, 6, 8, 10]
DUAL_VAR_BOUNDARY_CASES = [3, 7, 9]

MESH_STRESS_CASES = [
    {
        "id": f"case2_mesh_{n_elem}",
        "params": {
            "CASE_ID": 2,
            "F": -5000.0,
            "M_e": 0.0,
            "q": 0.0,
            "a": 0.5,
            "n_elem": n_elem,
        },
    }
    for n_elem in (2, 5, 10, 20, 50, 100, 200, 500)
]


def _as_jsonable(value):
    if isinstance(value, dict):
        return {key: _as_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_as_jsonable(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def make_params(overrides):
    params = config.PARAMS.copy()
    params.update(overrides)

    b = params.get("b", config.PARAMS["b"])
    h = params.get("h", config.PARAMS["h"])
    params["I"] = (b * h**3) / 12.0
    params["A"] = b * h
    params["BC_TYPE"] = "cantilever" if params["CASE_ID"] <= 4 else "simply_supported"
    return params


def dataframe_to_arrays(df):
    return {
        column: np.asarray(df[column].to_numpy(dtype=float), dtype=np.float64)
        for column in NUMERIC_COLUMNS
    }


def collect_solver_snapshot(params):
    case_id = params["CASE_ID"]
    fem = dataframe_to_arrays(get_fem_result(params))
    rk = dataframe_to_arrays(get_rk_result(params))

    solver = AnalysisFunc()
    method_name = f"generate_situation_{case_id}_data"
    solver_method = getattr(solver, method_name)
    ana_df = solver_method(
        E=params["E"],
        I=params["I"],
        l=params["L"],
        q=params.get("q", 0.0),
        F=params.get("F", 0.0),
        M_e=params.get("M_e", 0.0),
        a=params.get("a", 0.0),
        num_points=params["sample_n"],
    )
    ana = dataframe_to_arrays(ana_df)

    return {"fem": fem, "rk": rk, "ana": ana}


def collect_mesh_stress_snapshot(params):
    fem = dataframe_to_arrays(get_fem_result(params))
    max_abs_w = float(np.max(np.abs(fem["w"])))
    return {"fem": fem, "max_abs_w": max_abs_w}


def _make_analyzer(case_id, tolerance_percent=5.0):
    analyzer = AdaptiveFemAnalyzer(tolerance_percent=tolerance_percent)
    analyzer.p["CASE_ID"] = case_id
    analyzer.p["BC_TYPE"] = "cantilever" if case_id <= 4 else "simply_supported"
    return analyzer


def collect_single_var_boundary_snapshot(case_id, tolerance_percent=5.0):
    analyzer = _make_analyzer(case_id, tolerance_percent=tolerance_percent)
    conf = analyzer.case_config[case_id]
    load_var_name = conf["load"]

    sign = 1.0 if analyzer.p.get(load_var_name, -1.0) >= 0 else -1.0
    max_load = sign * 1.0

    p_test = analyzer.p.copy()
    for key in ("F", "M_e", "q"):
        p_test[key] = 0.0

    while True:
        p_test[load_var_name] = max_load
        p_test["CASE_ID"] = case_id
        err = analyzer._get_max_relative_error(p_test)

        if err >= analyzer.x0 * 1.2 or err == 999.0:
            break
        if abs(max_load) > 1e7:
            break
        max_load *= 1.5

    load_arr = np.linspace(0.0, max_load, 100)
    err_arr = []
    for val in load_arr:
        if val == 0:
            err_arr.append(0.0)
            continue
        p_test[load_var_name] = float(val)
        err = analyzer._get_max_relative_error(p_test)
        err_arr.append(float(err) if err != 999.0 else np.nan)

    exact_crit = float(
        analyzer.find_critical_load_fast(
            case_id,
            load_var_name,
            fixed_a=analyzer.p.get("L", 0.3) / 2.0,
        )
    )

    return {
        "load_var": load_var_name,
        "tolerance_percent": np.array([analyzer.x0], dtype=np.float64),
        "max_load": np.array([max_load], dtype=np.float64),
        "load_arr": np.asarray(load_arr, dtype=np.float64),
        "err_arr": np.asarray(err_arr, dtype=np.float64),
        "exact_crit": np.array([exact_crit], dtype=np.float64),
    }


def collect_dual_var_boundary_snapshot(
    case_id,
    max_anchor_points=90,
    chunk_size=5,
    workers="auto",
    start_method="forkserver",
):
    load_var, a_values = _case_metadata(case_id, max_anchor_points)
    tasks = _build_tasks([case_id], {case_id: a_values}, chunk_size)

    cpu_ids = _available_cpus()
    resolved_workers = _parse_workers(workers, cpu_ids)
    resolved_start_method = _choose_start_method(start_method)

    chunk_results = _run_parallel_tasks(
        tasks,
        resolved_workers,
        resolved_start_method,
        False,
        cpu_ids,
        verbose=False,
    )
    rows = _group_rows(chunk_results)[case_id]

    return {
        "load_var": load_var,
        "max_anchor_points": np.array([max_anchor_points], dtype=np.float64),
        "chunk_size": np.array([chunk_size], dtype=np.float64),
        "a_values": np.asarray([row[1] for row in rows], dtype=np.float64),
        "critical_loads": np.asarray([row[2] for row in rows], dtype=np.float64),
    }


def build_manifest():
    solver_cases = [
        {"id": case["id"], "params": make_params(case["p"])}
        for case in BATCH_TEST_CASES
    ]

    mesh_cases = [
        {"id": case["id"], "params": make_params(case["params"])}
        for case in MESH_STRESS_CASES
    ]

    return {
        "format_version": 1,
        "numeric_columns": list(NUMERIC_COLUMNS),
        "solver_cases": [
            {"id": case["id"], "params": _as_jsonable(case["params"])}
            for case in solver_cases
        ],
        "mesh_stress_cases": [
            {"id": case["id"], "params": _as_jsonable(case["params"])}
            for case in mesh_cases
        ],
        "single_var_boundary_cases": SINGLE_VAR_BOUNDARY_CASES,
        "dual_var_boundary_cases": DUAL_VAR_BOUNDARY_CASES,
        "dual_var_boundary_config": {
            "max_anchor_points": 90,
            "chunk_size": 5,
            "workers": "auto",
            "start_method": "forkserver",
        },
        "tolerances": {
            "solver_abs": 1e-8,
            "solver_rel": 1e-6,
            "mesh_abs": 1e-8,
            "mesh_rel": 1e-6,
            "boundary_abs": 1e-3,
            "boundary_rel": 1e-4,
        },
    }


def save_npz(path, arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    flat = {}

    def visit(prefix, value):
        if isinstance(value, dict):
            for key, item in value.items():
                next_prefix = f"{prefix}__{key}" if prefix else key
                visit(next_prefix, item)
            return
        flat[prefix] = np.asarray(value)

    visit("", arrays)
    np.savez_compressed(path, **flat)


def load_manifest():
    manifest_path = BASELINE_ROOT / "manifest.json"
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def load_npz(path):
    with np.load(path, allow_pickle=False) as data:
        root = {}
        for key in data.files:
            parts = key.split("__")
            current = root
            for part in parts[:-1]:
                current = current.setdefault(part, {})
            current[parts[-1]] = data[key]
        return root
