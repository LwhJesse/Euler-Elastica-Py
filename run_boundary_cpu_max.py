#!/usr/bin/env python3
"""
CPU-maximized scheduler for bi-variable FEM critical boundary scans.

This script preserves the numerical path from critical_boundary_fem.py and only
changes how load-position anchors are scheduled across CPU worker processes.
"""

import os

# Limit nested BLAS/OpenMP threading before importing NumPy/SciPy/OpenSees paths.
for _name in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/euler_elastica_matplotlib")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import argparse
import csv
import multiprocessing as mp
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from core import config

DEFAULT_CASES = [3, 7, 9]
TOLERANCE_PERCENT = config.BOUNDARY_TOLERANCE_PERCENT
OUTPUT_DIR = Path("results/boundary_analysis").resolve()

_WORKER_INDEX = None


def _available_cpus():
    if hasattr(os, "sched_getaffinity"):
        try:
            cpus = sorted(os.sched_getaffinity(0))
            if cpus:
                return cpus
        except OSError:
            pass

    count = os.cpu_count() or 1
    return list(range(count))


def _init_worker(pin_cores, counter, lock, cpu_ids):
    global _WORKER_INDEX

    worker_index = 0
    if counter is not None and lock is not None:
        with lock:
            worker_index = counter.value
            counter.value += 1

    _WORKER_INDEX = worker_index

    if not pin_cores:
        return

    if not hasattr(os, "sched_setaffinity"):
        print(
            "warning: CPU affinity requested but os.sched_setaffinity is unavailable",
            file=sys.stderr,
            flush=True,
        )
        return

    try:
        cpu = int(cpu_ids[worker_index % len(cpu_ids)])
        os.sched_setaffinity(0, {cpu})
    except Exception as exc:
        print(
            f"warning: failed to set CPU affinity for worker {worker_index}: {exc}",
            file=sys.stderr,
            flush=True,
        )


def _case_bc_type(case_id):
    return "cantilever" if case_id <= 4 else "simply_supported"


def _make_analyzer(case_id):
    from critical_boundary_fem import AdaptiveFemAnalyzer

    config.PARAMS["CASE_ID"] = case_id
    config.PARAMS["BC_TYPE"] = _case_bc_type(case_id)

    analyzer = AdaptiveFemAnalyzer(tolerance_percent=TOLERANCE_PERCENT)
    analyzer.p["CASE_ID"] = case_id
    analyzer.p["BC_TYPE"] = _case_bc_type(case_id)
    return analyzer


def _case_metadata(case_id, max_anchor_points):
    import numpy as np

    analyzer = _make_analyzer(case_id)
    conf = analyzer.case_config.get(case_id)
    if not conf:
        raise ValueError(f"Unknown CASE_ID = {case_id}")
    if not conf["two_var"]:
        raise ValueError(f"CASE_ID {case_id} is not a bi-variable boundary case")

    L = analyzer.p["L"]
    n_elem = analyzer.p.get("n_elem", 100)

    # Keep this logic identical to critical_boundary_fem.py's bi-variable branch.
    fem_nodes = np.linspace(0.0, L, n_elem + 1)
    valid_nodes = fem_nodes[(fem_nodes >= 0.05 * L) & (fem_nodes <= 0.95 * L)]

    if len(valid_nodes) > max_anchor_points:
        step = max(1, len(valid_nodes) // max_anchor_points)
        a_values = valid_nodes[::step]
    else:
        a_values = valid_nodes

    return conf["load"], [float(a) for a in a_values]


def _run_chunk(task):
    case_id = task["case_id"]
    chunk_id = task["chunk_id"]
    a_values = task["a_values"]

    started_at = time.time()
    analyzer = _make_analyzer(case_id)
    conf = analyzer.case_config[case_id]
    load_var = conf["load"]

    rows = []
    prev_crit = None
    for a in a_values:
        critical_load = analyzer.find_critical_load_fast(
            case_id,
            load_var,
            fixed_a=float(a),
            prev_crit=prev_crit,
        )
        critical_load = float(critical_load)
        rows.append((case_id, float(a), critical_load))
        prev_crit = critical_load

    finished_at = time.time()
    return {
        "case_id": case_id,
        "chunk_id": chunk_id,
        "rows": rows,
        "duration": finished_at - started_at,
        "started_at": started_at,
        "finished_at": finished_at,
        "pid": os.getpid(),
        "worker_index": _WORKER_INDEX,
    }


def _run_serial_verify_task(points_by_case):
    rows = []

    for case_id in sorted(points_by_case):
        analyzer = _make_analyzer(case_id)
        conf = analyzer.case_config[case_id]
        load_var = conf["load"]
        prev_crit = None

        for a in points_by_case[case_id]:
            critical_load = analyzer.find_critical_load_fast(
                case_id,
                load_var,
                fixed_a=float(a),
                prev_crit=prev_crit,
            )
            critical_load = float(critical_load)
            rows.append((case_id, float(a), critical_load))
            prev_crit = critical_load

    return rows


def _choose_start_method(requested):
    available = mp.get_all_start_methods()
    if requested in available:
        return requested

    fallback = "spawn" if "spawn" in available else available[0]
    print(
        f"warning: start method '{requested}' is unavailable; using '{fallback}'",
        file=sys.stderr,
    )
    return fallback


def _parse_workers(value, cpu_ids):
    if value == "auto":
        return max(1, len(cpu_ids))

    try:
        workers = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--workers must be 'auto' or an integer") from exc

    if workers < 1:
        raise argparse.ArgumentTypeError("--workers must be at least 1")
    return workers


def _executor(max_workers, start_method, pin_cores, cpu_ids):
    ctx = mp.get_context(start_method)
    counter = ctx.Value("i", 0)
    lock = ctx.Lock()
    return ProcessPoolExecutor(
        max_workers=max_workers,
        mp_context=ctx,
        initializer=_init_worker,
        initargs=(pin_cores, counter, lock, cpu_ids),
    )


def _build_tasks(cases, a_values_by_case, chunk_size):
    per_case_chunks = {}
    for case_id in cases:
        chunks = []
        a_values = a_values_by_case[case_id]
        for start in range(0, len(a_values), chunk_size):
            chunks.append(a_values[start : start + chunk_size])
        per_case_chunks[case_id] = chunks

    tasks = []
    max_chunks = max((len(chunks) for chunks in per_case_chunks.values()), default=0)
    for chunk_index in range(max_chunks):
        for case_id in cases:
            chunks = per_case_chunks[case_id]
            if chunk_index >= len(chunks):
                continue
            tasks.append(
                {
                    "case_id": case_id,
                    "chunk_id": f"case{case_id}_chunk{chunk_index:04d}",
                    "a_values": chunks[chunk_index],
                }
            )
    return tasks


def _run_parallel_tasks(tasks, workers, start_method, pin_cores, cpu_ids, verbose=True):
    results = []
    if not tasks:
        return results

    with _executor(workers, start_method, pin_cores, cpu_ids) as executor:
        futures = [executor.submit(_run_chunk, task) for task in tasks]
        total = len(futures)

        for done_count, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            results.append(result)

            if verbose:
                print(
                    f"[{done_count:>3}/{total}] {result['chunk_id']} "
                    f"finished in {result['duration']:.2f}s "
                    f"({len(result['rows'])} points, pid={result['pid']})",
                    flush=True,
                )

    return results


def _select_verify_points(a_values):
    if not a_values:
        return []

    indices = [0, len(a_values) // 2, len(a_values) - 1]
    selected = []
    seen = set()
    for index in indices:
        if index in seen:
            continue
        seen.add(index)
        selected.append(a_values[index])
    return selected


def _verify_small(cases, a_values_by_case, workers, start_method, pin_cores, cpu_ids):
    points_by_case = {
        case_id: _select_verify_points(a_values_by_case[case_id])
        for case_id in cases
    }
    total_points = sum(len(points) for points in points_by_case.values())
    if total_points == 0:
        print("verify-small: no points selected; skipping verification")
        return

    print(f"verify-small: computing {total_points} points serially")
    with _executor(1, start_method, False, cpu_ids) as executor:
        serial_rows = executor.submit(_run_serial_verify_task, points_by_case).result()

    verify_tasks = [
        {
            "case_id": case_id,
            "chunk_id": f"verify_case{case_id}",
            "a_values": points,
        }
        for case_id, points in points_by_case.items()
        if points
    ]
    verify_workers = min(max(1, workers), max(1, len(verify_tasks)))

    print(f"verify-small: computing {total_points} points through parallel chunk path")
    parallel_results = _run_parallel_tasks(
        verify_tasks,
        verify_workers,
        start_method,
        pin_cores,
        cpu_ids,
        verbose=False,
    )
    parallel_rows = [
        row
        for result in parallel_results
        for row in result["rows"]
    ]

    serial_map = {(case_id, a): critical for case_id, a, critical in serial_rows}
    parallel_map = {(case_id, a): critical for case_id, a, critical in parallel_rows}

    failures = 0
    for case_id in cases:
        for a in points_by_case[case_id]:
            key = (case_id, a)
            serial_crit = serial_map.get(key)
            parallel_crit = parallel_map.get(key)
            if serial_crit is None or parallel_crit is None:
                print(f"warning: verify-small missing result for case {case_id}, a={a:.12g}")
                failures += 1
                continue

            abs_err = abs(serial_crit - parallel_crit)
            rel_err = abs_err / max(abs(serial_crit), 1e-12)
            passed = abs_err <= 1e-3 or rel_err <= 1e-4

            if not passed:
                print(
                    "warning: verify-small mismatch "
                    f"case={case_id}, a={a:.12g}, "
                    f"serial={serial_crit:.12g}, parallel={parallel_crit:.12g}, "
                    f"abs_err={abs_err:.6g}, rel_err={rel_err:.6g}",
                    flush=True,
                )
                failures += 1

    if failures:
        print(f"verify-small: completed with {failures} warning(s)")
    else:
        print("verify-small: passed")


def _group_rows(chunk_results):
    rows_by_case = defaultdict(list)
    for result in chunk_results:
        rows_by_case[result["case_id"]].extend(result["rows"])

    for case_id in rows_by_case:
        rows_by_case[case_id].sort(key=lambda row: row[1])

    return rows_by_case


def _write_case_outputs(case_id, rows, load_var, suffix):
    import matplotlib.pyplot as plt

    if not rows:
        raise RuntimeError(f"No boundary rows produced for case {case_id}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    suffix_part = f"_{suffix}" if suffix else ""

    csv_path = OUTPUT_DIR / f"pure_critical_boundary_fem_case{case_id}{suffix_part}.csv"
    jpg_path = OUTPUT_DIR / f"pure_critical_boundary_fem_case{case_id}{suffix_part}.jpg"

    with csv_path.open("w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["case_id", "a", "critical_load"])
        for row_case_id, a, critical_load in rows:
            writer.writerow([row_case_id, f"{a:.12g}", f"{critical_load:.12g}"])

    a_values = [row[1] for row in rows]
    critical_loads = [row[2] for row in rows]

    plt.figure(figsize=(10, 6))
    plt.plot(
        a_values,
        critical_loads,
        "r-",
        lw=3,
        label=f"Critical Boundary ({TOLERANCE_PERCENT:g}% {config.ERROR_METRIC.get('label', 'Metric Error')})",
    )

    y_extreme = (
        min(critical_loads) * 1.2
        if min(critical_loads) < 0
        else max(critical_loads) * 1.2
    )
    plt.fill_between(
        a_values,
        critical_loads,
        y_extreme,
        color="#ffcccc",
        alpha=0.5,
        label="Nonlinear Zone (Danger)",
    )
    plt.fill_between(
        a_values,
        0,
        critical_loads,
        color="#e6f2e6",
        alpha=0.5,
        label="Linear Zone (Safe)",
    )

    plt.title(
        f"Case {case_id}: FEM Critical {load_var} vs Load Position $a$ "
        f"({config.ERROR_METRIC.get('label', 'Metric Error')})"
    )
    plt.xlabel("Load Position $a$ (m)")
    plt.ylabel(f"Critical Load {load_var}")
    plt.grid(True, linestyle="--")
    plt.legend()

    plt.savefig(jpg_path, dpi=300)
    plt.close()

    return csv_path, jpg_path


def _print_summary(
    total_elapsed,
    cases,
    rows_by_case,
    chunk_results,
    workers,
    chunk_size,
    output_paths,
):
    print("\nCPU boundary scan summary")
    print(f"Total elapsed: {total_elapsed:.2f}s ({total_elapsed / 60.0:.2f} min)")
    print(f"Workers: {workers}")
    print(f"Chunk size: {chunk_size}")
    print(f"Total chunks: {len(chunk_results)}")

    if chunk_results:
        avg_chunk = sum(result["duration"] for result in chunk_results) / len(chunk_results)
        print(f"Average chunk elapsed: {avg_chunk:.2f}s")
    else:
        print("Average chunk elapsed: n/a")

    for case_id in cases:
        case_chunks = [
            result for result in chunk_results
            if result["case_id"] == case_id
        ]
        if case_chunks:
            case_elapsed = (
                max(result["finished_at"] for result in case_chunks)
                - min(result["started_at"] for result in case_chunks)
            )
            chunk_cpu = sum(result["duration"] for result in case_chunks)
            print(
                f"Case {case_id}: {case_elapsed:.2f}s wall, "
                f"{chunk_cpu:.2f}s summed chunk time, "
                f"{len(rows_by_case[case_id])} points"
            )
        else:
            print(f"Case {case_id}: no chunks completed")

    print("Output files:")
    for csv_path, jpg_path in output_paths:
        print(f"  {csv_path}")
        print(f"  {jpg_path}")


def _parse_args(argv):
    parser = argparse.ArgumentParser(
        description="CPU-maximized dynamic chunk scheduler for CASE_ID 3/7/9 FEM boundaries."
    )
    parser.add_argument("--cases", nargs="+", type=int, default=DEFAULT_CASES)
    parser.add_argument("--max-anchor-points", type=int, default=90)
    parser.add_argument("--chunk-size", type=int, default=5)
    parser.add_argument("--workers", default="auto")

    pin_group = parser.add_mutually_exclusive_group()
    pin_group.add_argument("--pin-cores", dest="pin_cores", action="store_true")
    pin_group.add_argument("--no-pin-cores", dest="pin_cores", action="store_false")
    parser.set_defaults(pin_cores=False)

    parser.add_argument("--suffix", default="cpu_max")
    parser.add_argument(
        "--start-method",
        choices=("fork", "forkserver", "spawn"),
        default="forkserver",
    )
    parser.add_argument("--verify-small", action="store_true")

    args = parser.parse_args(argv)

    if args.max_anchor_points < 1:
        parser.error("--max-anchor-points must be at least 1")
    if args.chunk_size < 1:
        parser.error("--chunk-size must be at least 1")

    return args


def main(argv=None):
    script_started_at = time.time()
    args = _parse_args(argv)

    cases = list(dict.fromkeys(args.cases))
    cpu_ids = _available_cpus()
    try:
        workers = _parse_workers(args.workers, cpu_ids)
    except argparse.ArgumentTypeError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    start_method = _choose_start_method(args.start_method)

    load_vars_by_case = {}
    a_values_by_case = {}
    for case_id in cases:
        load_var, a_values = _case_metadata(case_id, args.max_anchor_points)
        load_vars_by_case[case_id] = load_var
        a_values_by_case[case_id] = a_values

    tasks = _build_tasks(cases, a_values_by_case, args.chunk_size)

    print("CPU-maximized FEM boundary scan")
    print(f"Cases: {cases}")
    print(f"Workers: {workers} ({args.workers})")
    print(f"Start method: {start_method}")
    print(f"Pin cores: {args.pin_cores}")
    print(f"Chunk size: {args.chunk_size}")
    print(f"Total chunks: {len(tasks)}")
    for case_id in cases:
        print(
            f"Case {case_id}: {len(a_values_by_case[case_id])} anchor points, "
            f"load variable={load_vars_by_case[case_id]}"
        )

    if args.verify_small:
        _verify_small(
            cases,
            a_values_by_case,
            workers,
            start_method,
            args.pin_cores,
            cpu_ids,
        )

    chunk_results = _run_parallel_tasks(
        tasks,
        workers,
        start_method,
        args.pin_cores,
        cpu_ids,
        verbose=True,
    )
    rows_by_case = _group_rows(chunk_results)

    output_paths = []
    for case_id in cases:
        rows = rows_by_case[case_id]
        expected_count = len(a_values_by_case[case_id])
        if len(rows) != expected_count:
            print(
                f"warning: case {case_id} produced {len(rows)} rows, "
                f"expected {expected_count}",
                file=sys.stderr,
            )
        output_paths.append(
            _write_case_outputs(
                case_id,
                rows,
                load_vars_by_case[case_id],
                args.suffix,
            )
        )

    total_elapsed = time.time() - script_started_at
    _print_summary(
        total_elapsed,
        cases,
        rows_by_case,
        chunk_results,
        workers,
        args.chunk_size,
        output_paths,
    )
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    raise SystemExit(main())
