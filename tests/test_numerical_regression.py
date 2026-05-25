import math
import unittest
from pathlib import Path

import numpy as np

from tests.numerical_regression_cases import (
    BASELINE_ROOT,
    NUMERIC_COLUMNS,
    collect_dual_var_boundary_snapshot,
    collect_mesh_stress_snapshot,
    collect_single_var_boundary_snapshot,
    collect_solver_snapshot,
    load_manifest,
    load_npz,
    make_params,
)


def assert_close_or_relative(test_case, actual, expected, abs_tol, rel_tol, label):
    actual_arr = np.atleast_1d(np.asarray(actual, dtype=np.float64))
    expected_arr = np.atleast_1d(np.asarray(expected, dtype=np.float64))

    test_case.assertEqual(
        actual_arr.shape,
        expected_arr.shape,
        f"{label}: shape mismatch {actual_arr.shape} != {expected_arr.shape}",
    )

    diff = np.abs(actual_arr - expected_arr)
    rel = diff / np.maximum(np.abs(expected_arr), 1e-12)
    ok = np.isfinite(diff) & ((diff <= abs_tol) | (rel <= rel_tol))

    nan_match = np.isnan(actual_arr) & np.isnan(expected_arr)
    ok = ok | nan_match

    if np.all(ok):
        return

    bad_indices = np.argwhere(~ok)
    first_bad = tuple(int(index) for index in bad_indices[0])
    test_case.fail(
        f"{label}: mismatch at index {first_bad}, "
        f"actual={actual_arr[first_bad]!r}, expected={expected_arr[first_bad]!r}, "
        f"abs_diff={diff[first_bad]!r}, rel_diff={rel[first_bad]!r}"
    )


class NumericalRegressionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = load_manifest()
        cls.tolerances = cls.manifest["tolerances"]

    def test_baselines_exist(self):
        self.assertTrue((BASELINE_ROOT / "manifest.json").exists())
        self.assertTrue((BASELINE_ROOT / "solver_cases").exists())
        self.assertTrue((BASELINE_ROOT / "mesh_stress_cases").exists())
        self.assertTrue((BASELINE_ROOT / "single_var_boundaries").exists())
        self.assertTrue((BASELINE_ROOT / "dual_var_boundaries").exists())

    def test_solver_regressions(self):
        abs_tol = self.tolerances["solver_abs"]
        rel_tol = self.tolerances["solver_rel"]

        for case in self.manifest["solver_cases"]:
            case_id = case["id"]
            with self.subTest(case=case_id):
                expected = load_npz(BASELINE_ROOT / "solver_cases" / f"{case_id}.npz")
                actual = collect_solver_snapshot(make_params(case["params"]))

                for solver_name in ("fem", "rk", "ana"):
                    for column in NUMERIC_COLUMNS:
                        label = f"{case_id}:{solver_name}:{column}"
                        assert_close_or_relative(
                            self,
                            actual[solver_name][column],
                            expected[solver_name][column],
                            abs_tol,
                            rel_tol,
                            label,
                        )

    def test_mesh_stress_regressions(self):
        abs_tol = self.tolerances["mesh_abs"]
        rel_tol = self.tolerances["mesh_rel"]

        for case in self.manifest["mesh_stress_cases"]:
            case_id = case["id"]
            with self.subTest(case=case_id):
                expected = load_npz(BASELINE_ROOT / "mesh_stress_cases" / f"{case_id}.npz")
                actual = collect_mesh_stress_snapshot(make_params(case["params"]))

                for column in NUMERIC_COLUMNS:
                    assert_close_or_relative(
                        self,
                        actual["fem"][column],
                        expected["fem"][column],
                        abs_tol,
                        rel_tol,
                        f"{case_id}:fem:{column}",
                    )

                assert_close_or_relative(
                    self,
                    np.array([actual["max_abs_w"]], dtype=np.float64),
                    expected["max_abs_w"],
                    abs_tol,
                    rel_tol,
                    f"{case_id}:max_abs_w",
                )

    def test_single_var_boundary_regressions(self):
        abs_tol = self.tolerances["boundary_abs"]
        rel_tol = self.tolerances["boundary_rel"]

        for case_id in self.manifest["single_var_boundary_cases"]:
            with self.subTest(case=case_id):
                expected = load_npz(BASELINE_ROOT / "single_var_boundaries" / f"case{case_id}.npz")
                actual = collect_single_var_boundary_snapshot(case_id)

                for key in ("max_load", "load_arr", "err_arr", "exact_crit"):
                    assert_close_or_relative(
                        self,
                        actual[key],
                        expected[key],
                        abs_tol,
                        rel_tol,
                        f"single_var_case{case_id}:{key}",
                    )

                self.assertEqual(actual["load_var"], expected["load_var"].item())

    def test_dual_var_boundary_regressions(self):
        abs_tol = self.tolerances["boundary_abs"]
        rel_tol = self.tolerances["boundary_rel"]
        config = self.manifest["dual_var_boundary_config"]

        for case_id in self.manifest["dual_var_boundary_cases"]:
            with self.subTest(case=case_id):
                expected = load_npz(BASELINE_ROOT / "dual_var_boundaries" / f"case{case_id}.npz")
                actual = collect_dual_var_boundary_snapshot(case_id, **config)

                for key in ("max_anchor_points", "chunk_size", "a_values", "critical_loads"):
                    assert_close_or_relative(
                        self,
                        actual[key],
                        expected[key],
                        abs_tol,
                        rel_tol,
                        f"dual_var_case{case_id}:{key}",
                    )

                self.assertEqual(actual["load_var"], expected["load_var"].item())


if __name__ == "__main__":
    unittest.main(verbosity=2)
