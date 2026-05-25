import json
from pathlib import Path

from tests.numerical_regression_cases import (
    BASELINE_ROOT,
    build_manifest,
    collect_dual_var_boundary_snapshot,
    collect_mesh_stress_snapshot,
    collect_single_var_boundary_snapshot,
    collect_solver_snapshot,
    make_params,
    save_npz,
)


def main():
    manifest = build_manifest()
    BASELINE_ROOT.mkdir(parents=True, exist_ok=True)

    for case in manifest["solver_cases"]:
        snapshot = collect_solver_snapshot(make_params(case["params"]))
        save_npz(BASELINE_ROOT / "solver_cases" / f"{case['id']}.npz", snapshot)
        print(f"generated solver baseline: {case['id']}", flush=True)

    for case in manifest["mesh_stress_cases"]:
        snapshot = collect_mesh_stress_snapshot(make_params(case["params"]))
        save_npz(BASELINE_ROOT / "mesh_stress_cases" / f"{case['id']}.npz", snapshot)
        print(f"generated mesh-stress baseline: {case['id']}", flush=True)

    for case_id in manifest["single_var_boundary_cases"]:
        snapshot = collect_single_var_boundary_snapshot(case_id)
        save_npz(BASELINE_ROOT / "single_var_boundaries" / f"case{case_id}.npz", snapshot)
        print(f"generated single-var boundary baseline: case{case_id}", flush=True)

    boundary_config = manifest["dual_var_boundary_config"]
    for case_id in manifest["dual_var_boundary_cases"]:
        snapshot = collect_dual_var_boundary_snapshot(case_id, **boundary_config)
        save_npz(BASELINE_ROOT / "dual_var_boundaries" / f"case{case_id}.npz", snapshot)
        print(f"generated dual-var boundary baseline: case{case_id}", flush=True)

    (BASELINE_ROOT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(f"wrote manifest: {BASELINE_ROOT / 'manifest.json'}", flush=True)


if __name__ == "__main__":
    main()
