# compare_data.py
# 一键对比两组序列型工程数据（如梁的 w/theta/M/V/N），自动识别结构，插值对齐，计算指标，并输出单张对比图。
# 入口：one_click_compare(left, right, ...)
# - left/right 支持：pandas.DataFrame / dict(含 final_df/nodes_df/sections_df/curve_df) / 路径(csv/json/parquet)
# - 自动识别坐标列：优先 x/s；自动识别数值列：优先 w/theta/M/V/N/u
# - 对齐：公共区间上用等距步长 dx 线性插值
# - 指标：RMSE/MAE/Max|err|/MeanErr/MAPE/NRMSE(range)/R2/Pearson-r/L1积/ L2积/斜率RMSE/MAE
# - 绘图：单张 PNG。每个字段一行：左图叠加(FEM vs RK)，右图残差(resid)。
# - 返回：dict，含 summary_df/overall_df/align_df/structure_df

import os
import json
import math
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

EPS = 1e-12

# ----------------------------
# 加载与检测
# ----------------------------

def _load_any(obj: Any) -> pd.DataFrame:
    # 支持 DataFrame / dict / 文件路径
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    if isinstance(obj, dict):
        # 常见键的优先级：final_df > curve_df > nodes_df > sections_df > reactions_df
        for k in ["final_df", "curve_df", "nodes_df", "sections_df", "reactions_df", "data", "df"]:
            if k in obj and isinstance(obj[k], pd.DataFrame):
                return obj[k].copy()
        # 尝试把 dict 当表
        try:
            return pd.DataFrame(obj).copy()
        except Exception as _:
            raise ValueError("dict 中未找到 DataFrame（final_df/curve_df/nodes_df/sections_df/reactions_df/data/df）")
    if isinstance(obj, str):
        path = obj
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        ext = os.path.splitext(path)[1].lower()
        if ext in [".csv", ".tsv"]:
            sep = "," if ext == ".csv" else "\t"
            return pd.read_csv(path, sep=sep)
        if ext in [".json"]:
            try:
                with open(path, "r") as f:
                    payload = json.load(f)
                if isinstance(payload, list):
                    return pd.DataFrame(payload)
                if isinstance(payload, dict):
                    for k in ["final_df", "curve_df", "nodes_df", "sections_df", "reactions_df", "data", "df"]:
                        if k in payload and isinstance(payload[k], list):
                            return pd.DataFrame(payload[k])
                    return pd.DataFrame(payload)
            except Exception:
                return pd.read_json(path, lines=False)
        if ext in [".parquet", ".pq"]:
            return pd.read_parquet(path)
        raise ValueError("不支持的文件类型: %s" % ext)
    raise TypeError("输入必须是 DataFrame / dict / 路径字符串")

def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return df

def _detect_x_col(df: pd.DataFrame, prefer: Optional[str] = None) -> str:
    cols = set(df.columns)
    if prefer and prefer.lower() in cols:
        return prefer.lower()
    for c in ["x", "s", "dist", "position"]:
        if c in cols:
            return c
    num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if not num_cols:
        raise ValueError("找不到数值列作为坐标")
    # 选单调性较好的列
    best = None
    best_score = -1.0
    for c in num_cols:
        v = df[c].to_numpy()
        if len(v) < 2:
            continue
        d = np.diff(v)
        score = np.mean(d >= -1e-15)
        if score > best_score:
            best_score = score
            best = c
    return best or num_cols[0]

def _detect_value_cols(df: pd.DataFrame, x_col: str) -> List[str]:
    priority = ["w", "u", "theta", "m", "v", "n", "rz"]
    cols = [c for c in priority if c in df.columns]
    if not cols:
        cols = [c for c in df.columns if c != x_col and pd.api.types.is_numeric_dtype(df[c])]
    # 去重保序
    seen = set(); out = []
    for c in cols:
        if c not in seen:
            seen.add(c); out.append(c)
    return out

def _coerce_numeric(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    df = df.copy()
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

# ----------------------------
# 对齐与插值
# ----------------------------

def _overlap_uniform_grid(xl: np.ndarray, xr: np.ndarray, dx: float) -> np.ndarray:
    xmin = max(np.nanmin(xl), np.nanmin(xr))
    xmax = min(np.nanmax(xl), np.nanmax(xr))
    if xmax - xmin < 10 * EPS:
        raise ValueError("两组数据的 x 区间无重叠")
    n = max(2, int(round((xmax - xmin) / dx)) + 1)
    return np.linspace(xmin, xmax, n)

def _interp_to(xs: np.ndarray, xsrc: np.ndarray, ysrc: np.ndarray) -> np.ndarray:
    # 线性插值，范围外 NaN
    order = np.argsort(xsrc)
    x1 = xsrc[order]; y1 = ysrc[order]
    mask = np.isfinite(x1) & np.isfinite(y1)
    x1 = x1[mask]; y1 = y1[mask]
    if len(x1) == 0:
        return np.full_like(xs, np.nan, dtype=float)
    y = np.interp(xs, x1, y1)
    y[(xs < x1.min()) | (xs > x1.max())] = np.nan
    return y

def _derivative_cdiff(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    dydx = np.empty_like(y)
    if x.size < 3:
        slope = (y[-1] - y[0]) / max(EPS, (x[-1] - x[0]))
        dydx.fill(slope); return dydx
    dx = np.diff(x)
    dydx[1:-1] = (y[2:] - y[:-2]) / (x[2:] - x[:-2])
    dydx[0] = (y[1] - y[0]) / max(EPS, dx[0])
    dydx[-1] = (y[-1] - y[-2]) / max(EPS, dx[-1])
    return dydx

# ----------------------------
# 指标
# ----------------------------

def _safe_mape(y_true: np.ndarray, y_pred: np.ndarray, eps_abs: float = 1e-9) -> float:
    m = np.abs(y_true) > eps_abs
    if not np.any(m):
        return float('nan')
    return float(np.mean(np.abs((y_pred[m] - y_true[m]) / y_true[m])))

def _metrics(x: np.ndarray, y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    m = np.isfinite(y_true) & np.isfinite(y_pred)
    if not np.any(m):
        return {k: float('nan') for k in [
            "count","rmse","mae","max_abs_err","mean_err","mape",
            "nrmse_range","r2","pearson_r","area_abs_resid","area_sq_resid",
            "area_rel_L1","slope_rmse","slope_mae"
        ]}
    xt = x[m]; yt = y_true[m]; yp = y_pred[m]
    resid = yp - yt
    sse = float(np.dot(resid, resid))
    mse = sse / yt.size
    rmse = math.sqrt(mse)
    mae = float(np.mean(np.abs(resid)))
    max_abs_err = float(np.max(np.abs(resid)))
    mean_err = float(np.mean(resid))
    ybar = float(np.mean(yt))
    sst = float(np.dot(yt - ybar, yt - ybar))
    r2 = float(1.0 - sse / sst) if sst > EPS else float('nan')
    pearson_r = float(np.corrcoef(yt, yp)[0, 1]) if (np.std(yt) * np.std(yp) > EPS) else float('nan')
    yrng = float(np.max(yt) - np.min(yt))
    nrmse_range = float(rmse / yrng) if yrng > EPS else float('nan')
    area_abs_resid = float(np.trapz(np.abs(resid), xt))
    area_sq_resid = float(np.trapz(resid * resid, xt))
    area_true_L1 = float(np.trapz(np.abs(yt), xt))
    area_rel_L1 = float(area_abs_resid / area_true_L1) if area_true_L1 > EPS else float('nan')
    dy_true = _derivative_cdiff(xt, yt)
    dy_pred = _derivative_cdiff(xt, yp)
    dy_res = dy_pred - dy_true
    slope_rmse = float(math.sqrt(np.mean(dy_res * dy_res)))
    slope_mae = float(np.mean(np.abs(dy_res)))
    mape = _safe_mape(yt, yp)
    return {
        "count": int(yt.size),
        "rmse": rmse,
        "mae": mae,
        "max_abs_err": max_abs_err,
        "mean_err": mean_err,
        "mape": mape,
        "nrmse_range": nrmse_range,
        "r2": r2,
        "pearson_r": pearson_r,
        "area_abs_resid": area_abs_resid,
        "area_sq_resid": area_sq_resid,
        "area_rel_L1": area_rel_L1,
        "slope_rmse": slope_rmse,
        "slope_mae": slope_mae,
    }

# ----------------------------
# 主函数：一键对比 + 单张图
# ----------------------------

def one_click_compare(
    left: Union[pd.DataFrame, dict, str],
    right: Union[pd.DataFrame, dict, str],
    *,
    x_left: Optional[str] = None,
    x_right: Optional[str] = None,
    fields: Optional[List[str]] = None,   # 指定要比对的列；None 则自动挑选共有列
    dx: float = 0.0005,                   # 重采样步长
    out_png: str = "compare.png",
    title: str = "Data Comparison",
    left_label: str = "RK",
    right_label: str = "FEM",
    save_prefix: Optional[str] = None     # 若给定，将另存 CSV（summary/overall/align/structure）
) -> Dict[str, pd.DataFrame]:
    """
    一键对比两组数据，输出单张图和各类指标。
    返回 dict：{summary_df, overall_df, align_df, structure_df}
    """
    # 1) 加载与规范化
    dfL = _normalize_columns(_load_any(left))
    dfR = _normalize_columns(_load_any(right))

    # 2) 坐标列识别
    xL = _detect_x_col(dfL, x_left)
    xR = _detect_x_col(dfR, x_right)

    # 3) 值列识别
    colsL = _detect_value_cols(dfL, xL)
    colsR = _detect_value_cols(dfR, xR)
    # 常见字段优先
    pref = ["w", "theta", "m", "v", "n", "u"]
    common = [c for c in pref if (c in colsL and c in colsR)]
    if fields:
        fields = [c.lower() for c in fields if (c.lower() in dfL.columns and c.lower() in dfR.columns)]
    else:
        fields = common if common else sorted(list(set(colsL).intersection(set(colsR))))
    if not fields:
        raise ValueError("两组数据没有可以比较的公共数值列")

    # 4) 数值化
    dfL = _coerce_numeric(dfL, [xL] + fields)
    dfR = _coerce_numeric(dfR, [xR] + fields)

    # 5) 公共等距网格
    xs = _overlap_uniform_grid(dfL[xL].to_numpy(float), dfR[xR].to_numpy(float), dx=dx)

    # 6) 插值、指标
    align_cols: Dict[str, np.ndarray] = {"x": xs}
    rows = []
    for c in fields:
        yL = _interp_to(xs, dfL[xL].to_numpy(float), dfL[c].to_numpy(float))
        yR = _interp_to(xs, dfR[xR].to_numpy(float), dfR[c].to_numpy(float))
        align_cols[f"{c}_left"] = yL
        align_cols[f"{c}_right"] = yR
        align_cols[f"{c}_diff"] = yR - yL
        met = _metrics(xs, yL, yR)
        row = {"field": c}; row.update(met)
        rows.append(row)

    align_df = pd.DataFrame(align_cols)
    summary_df = pd.DataFrame(rows).set_index("field").sort_index()

    overall = {
        "fields": ",".join(fields),
        "points": int(len(align_df)),
        "rmse_mean": float(np.nanmean(summary_df["rmse"])) if len(summary_df) else np.nan,
        "mae_mean": float(np.nanmean(summary_df["mae"])) if len(summary_df) else np.nan,
        "max_abs_max": float(np.nanmax(summary_df["max_abs_err"])) if len(summary_df) else np.nan,
        "nrmse_range_mean": float(np.nanmean(summary_df["nrmse_range"])) if len(summary_df) else np.nan,
        "r2_mean": float(np.nanmean(summary_df["r2"])) if len(summary_df) else np.nan,
        "pearson_r_mean": float(np.nanmean(summary_df["pearson_r"])) if len(summary_df) else np.nan,
        "mape_mean": float(np.nanmean(summary_df["mape"])) if len(summary_df) else np.nan,
        "slope_rmse_mean": float(np.nanmean(summary_df["slope_rmse"])) if len(summary_df) else np.nan,
    }
    overall_df = pd.DataFrame([overall])

    structure_df = pd.DataFrame([
        {"side": "left", "label": left_label, "x_col": xL, "value_cols": ",".join(fields),
         "n_rows": len(dfL), "dtypes": json.dumps({c: str(dfL[c].dtype) for c in dfL.columns})},
        {"side": "right", "label": right_label, "x_col": xR, "value_cols": ",".join(fields),
         "n_rows": len(dfR), "dtypes": json.dumps({c: str(dfR[c].dtype) for c in dfR.columns})},
    ])

    # 7) 单张图：每个字段一行，两列（左=叠加，右=残差）
    nrows = len(fields); ncols = 2
    h_each = 2.4  # 每行高度
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(10, max(3.5, h_each * nrows)), sharex=True)
    if nrows == 1:
        axes = np.array([axes])  # 统一二维索引

    for i, c in enumerate(fields):
        ax_curv = axes[i, 0]
        ax_res = axes[i, 1]
        ax_curv.plot(align_df["x"], align_df[f"{c}_left"], lw=2, label=left_label)
        ax_curv.plot(align_df["x"], align_df[f"{c}_right"], lw=2, ls="--", label=right_label)
        ax_curv.set_ylabel(c)
        ax_curv.grid(True, alpha=0.3)
        ax_curv.legend(loc="best", fontsize=9)

        resid = align_df[f"{c}_diff"]
        ax_res.plot(align_df["x"], resid, color="tab:red", lw=1.5)
        ax_res.axhline(0.0, color="k", lw=1)
        ax_res.set_ylabel(f"{c} resid")
        ax_res.grid(True, alpha=0.3)

        # 在叠加图角落放核心指标
        m = summary_df.loc[c].to_dict()
        box = (
            f"RMSE={m['rmse']:.6g}\nMAE={m['mae']:.6g}\n"
            f"NRMSE={m['nrmse_range']:.6g}\nR2={m['r2']:.6g}\n"
            f"r={m['pearson_r']:.6g}\nMAPE={m['mape']:.6g}"
        )
        ax_curv.text(0.02, 0.98, box, transform=ax_curv.transAxes, va="top",
                     bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), fontsize=8)

    axes[-1, 0].set_xlabel("x")
    axes[-1, 1].set_xlabel("x")
    fig.suptitle(title, y=0.995, fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

    # 8) 可选保存
    if save_prefix:
        summary_df.to_csv(f"{save_prefix}_summary.csv")
        overall_df.to_csv(f"{save_prefix}_overall.csv", index=False)
        align_df.to_csv(f"{save_prefix}_align.csv", index=False)
        structure_df.to_csv(f"{save_prefix}_structure.csv", index=False)

    return {
        "summary_df": summary_df,
        "overall_df": overall_df,
        "align_df": align_df,
        "structure_df": structure_df,
    }

# 可选：简单自测
if __name__ == "__main__":
    x = np.linspace(0, 1, 1001)
    rk = pd.DataFrame({"x": x, "w": -0.1*x, "theta": -0.2*x})
    fem = pd.DataFrame({"x": x, "w": -0.1*x + 0.002*np.sin(25*x), "theta": -0.2*x*1.01})
    res = one_click_compare(rk, fem, dx=0.001, out_png="compare_demo.png", title="Demo RK vs FEM")
    print(res["summary_df"])
    print(res["overall_df"])
