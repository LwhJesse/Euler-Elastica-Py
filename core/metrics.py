"""
Configurable error metrics for boundary-search data analysis.

All metrics compare already-generated solution DataFrames on a shared initial
arc-length coordinate ``s``. Solver code must not rely on row order matching.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np


def _integrate(y, x):
    if hasattr(np, "trapezoid"):
        return np.trapezoid(y, x)
    return np.trapz(y, x)


def _as_float(value, default):
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def _require_columns(df, required: Iterable[str], name: str):
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"{name} solution is missing required column(s): {missing}")


def _sorted_unique_xy(df, column):
    s = np.asarray(df["s"], dtype=float)
    y = np.asarray(df[column], dtype=float)
    valid = np.isfinite(s) & np.isfinite(y)
    if np.count_nonzero(valid) < 2:
        raise ValueError(f"column '{column}' does not contain enough finite samples")

    s = s[valid]
    y = y[valid]
    order = np.argsort(s)
    s = s[order]
    y = y[order]

    unique_s, inverse = np.unique(s, return_inverse=True)
    if unique_s.size < 2:
        raise ValueError("solution column 's' must contain at least two unique values")
    if unique_s.size == s.size:
        return unique_s, y

    y_accum = np.zeros_like(unique_s, dtype=float)
    counts = np.zeros_like(unique_s, dtype=float)
    np.add.at(y_accum, inverse, y)
    np.add.at(counts, inverse, 1.0)
    return unique_s, y_accum / counts


def align_solution_pair(df_ref, df_cmp, L, columns, grid_n=5000):
    """
    Interpolate two solution DataFrames to the same ``s`` grid.

    ``s`` is the only alignment key. Row index or row count are intentionally
    ignored because analytical, RK, and FEM paths sample the beam differently.
    """
    if not columns:
        raise ValueError("metric requires at least one comparison column")

    columns = list(columns)
    _require_columns(df_ref, ["s", *columns], "reference")
    _require_columns(df_cmp, ["s", *columns], "candidate")

    L = float(L)
    grid_n = int(grid_n)
    if L <= 0.0:
        raise ValueError(f"beam length L must be positive, got {L}")
    if grid_n < 2:
        raise ValueError(f"grid_n must be at least 2, got {grid_n}")

    s_grid = np.linspace(0.0, L, grid_n)
    ref_arrays = {}
    cmp_arrays = {}

    for column in columns:
        s_ref, y_ref = _sorted_unique_xy(df_ref, column)
        s_cmp, y_cmp = _sorted_unique_xy(df_cmp, column)
        ref_arrays[column] = np.interp(s_grid, s_ref, y_ref)
        cmp_arrays[column] = np.interp(s_grid, s_cmp, y_cmp)

    return s_grid, ref_arrays, cmp_arrays


def compute_error_metric(df_ref, df_cmp, L, metric_config):
    """Return the configured error metric as a percentage."""
    if not isinstance(metric_config, dict):
        raise TypeError("metric_config must be a dictionary")

    name = metric_config.get("name")
    if name == "industrial_centerline_shape":
        value = _industrial_centerline_shape(df_ref, df_cmp, L, metric_config)
    elif name == "max_relative":
        value = _max_relative(df_ref, df_cmp, L, metric_config)
    elif name == "l2_curve":
        value = _l2_curve(df_ref, df_cmp, L, metric_config)
    elif name == "l2_vector_curve":
        value = _l2_vector_curve(df_ref, df_cmp, L, metric_config)
    else:
        raise ValueError(f"unsupported error metric: {name!r}")

    if not np.isfinite(value):
        raise ValueError(f"metric {name!r} produced a non-finite value: {value}")
    return float(value)


def _columns(metric_config):
    cols = metric_config.get("columns")
    if not cols:
        raise ValueError("metric_config must define non-empty 'columns'")
    return list(cols)


def _epsilon(metric_config):
    return _as_float(metric_config.get("epsilon"), 1e-14)


def _grid_n(metric_config):
    return int(metric_config.get("grid_n", 5000))


def _stack(arrays, columns):
    return np.column_stack([arrays[column] for column in columns])


def _row_norm(values):
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        return np.abs(values)
    return np.linalg.norm(values, axis=1)


def _max_relative(df_ref, df_cmp, L, metric_config):
    columns = _columns(metric_config)
    eps = _epsilon(metric_config)
    s_grid, ref, cmp = align_solution_pair(df_ref, df_cmp, L, columns, _grid_n(metric_config))

    del s_grid
    ref_values = _stack(ref, columns)
    cmp_values = _stack(cmp, columns)
    diff_norm = _row_norm(cmp_values - ref_values)
    ref_norm = _row_norm(ref_values)

    normalization = metric_config.get("normalization", "reference_max_abs")
    if normalization != "reference_max_abs":
        raise ValueError(f"max_relative does not support normalization={normalization!r}")

    return 100.0 * np.max(diff_norm) / (np.max(ref_norm) + eps)


def _l2_curve(df_ref, df_cmp, L, metric_config):
    columns = _columns(metric_config)
    eps = _epsilon(metric_config)
    s_grid, ref, cmp = align_solution_pair(df_ref, df_cmp, L, columns, _grid_n(metric_config))

    diff_values = _stack(cmp, columns) - _stack(ref, columns)
    ref_values = _stack(ref, columns)
    diff_sq = np.sum(diff_values * diff_values, axis=1)
    ref_sq = np.sum(ref_values * ref_values, axis=1)

    normalization = metric_config.get("normalization", "reference_l2")
    if normalization != "reference_l2":
        raise ValueError(f"l2_curve does not support normalization={normalization!r}")

    return 100.0 * np.sqrt(_integrate(diff_sq, s_grid) / (_integrate(ref_sq, s_grid) + eps))


def _l2_vector_curve(df_ref, df_cmp, L, metric_config):
    return _l2_curve(df_ref, df_cmp, L, metric_config)


def _curvature(s_grid, arrays, metric_config):
    source = metric_config.get("curvature_source", "theta_gradient")
    if source == "theta_gradient":
        theta = np.asarray(arrays["theta"], dtype=float)
        return np.gradient(theta, s_grid, edge_order=1)
    if source == "moment_over_ei":
        if "M" not in arrays:
            raise ValueError("curvature_source='moment_over_ei' requires column 'M'")
        ei = metric_config.get("EI")
        if ei is None:
            raise ValueError("curvature_source='moment_over_ei' requires metric_config['EI']")
        return np.asarray(arrays["M"], dtype=float) / (float(ei) + _epsilon(metric_config))
    raise ValueError(f"unsupported curvature_source: {source!r}")


def _industrial_centerline_shape(df_ref, df_cmp, L, metric_config):
    eps = _epsilon(metric_config)
    grid_n = _grid_n(metric_config)
    columns = ["x", "w", "theta"]
    if metric_config.get("curvature_source") == "moment_over_ei":
        columns.append("M")
    s_grid, ref, cmp = align_solution_pair(df_ref, df_cmp, L, columns, grid_n)

    ref_r = np.column_stack([ref["x"], ref["w"]])
    cmp_r = np.column_stack([cmp["x"], cmp["w"]])
    r0 = np.column_stack([s_grid, np.zeros_like(s_grid)])

    ref_theta = np.asarray(ref["theta"], dtype=float)
    cmp_theta = np.asarray(cmp["theta"], dtype=float)
    ref_t = np.column_stack([np.cos(ref_theta), np.sin(ref_theta)])
    cmp_t = np.column_stack([np.cos(cmp_theta), np.sin(cmp_theta)])
    t0 = np.column_stack([np.ones_like(s_grid), np.zeros_like(s_grid)])

    ref_kappa = _curvature(s_grid, ref, metric_config)
    cmp_kappa = _curvature(s_grid, cmp, metric_config)

    er = np.sqrt(
        _integrate(np.sum((cmp_r - ref_r) ** 2, axis=1), s_grid)
        / (_integrate(np.sum((ref_r - r0) ** 2, axis=1), s_grid) + eps)
    )
    et = np.sqrt(
        _integrate(np.sum((cmp_t - ref_t) ** 2, axis=1), s_grid)
        / (_integrate(np.sum((ref_t - t0) ** 2, axis=1), s_grid) + eps)
    )
    ek = np.sqrt(
        _integrate((cmp_kappa - ref_kappa) ** 2, s_grid)
        / (_integrate(ref_kappa ** 2, s_grid) + eps)
    )

    weights = metric_config.get("weights", {})
    ar = _as_float(weights.get("position"), 1.0)
    at = _as_float(weights.get("tangent"), 0.5)
    ak = _as_float(weights.get("curvature"), 0.25)
    weight_sum = ar + at + ak
    if weight_sum <= 0.0:
        raise ValueError("industrial_centerline_shape weights must sum to a positive value")

    e_global = np.sqrt((ar * er * er + at * et * et + ak * ek * ek) / weight_sum)

    er_inf = np.max(_row_norm(cmp_r - ref_r)) / (np.max(_row_norm(ref_r - r0)) + eps)
    et_inf = np.max(_row_norm(cmp_t - ref_t)) / (np.max(_row_norm(ref_t - t0)) + eps)
    ek_inf = np.max(np.abs(cmp_kappa - ref_kappa)) / (np.max(np.abs(ref_kappa)) + eps)

    local_weights = metric_config.get("local_weights", {})
    lt = _as_float(local_weights.get("tangent"), 0.75)
    lk = _as_float(local_weights.get("curvature"), 0.5)
    e_local = max(er_inf, lt * et_inf, lk * ek_inf)

    return 100.0 * max(e_global, e_local)
