import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import interp1d
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import pandas as pd
# 引入你的求解器
import opensees_beam
import RK
import analysis_solve

# 参数设置
L = 1.0
E = 2e+11
I = (0.02*0.02**3)/12
A = 0.02*0.002

def compare_three_datasets(data_FEM, data_RK, data_Analysis, save_name='comparison_3sets.png'):

    # --- 1. 数据提取 (保持不变) ---
    def _resolve_df(obj):
        if isinstance(obj, dict):
            for k in ('final_df', 'nodes_df', 'sections_df'):
                if k in obj:
                    cand = obj[k]
                    if hasattr(cand, 'columns'):
                        return cand
            raise ValueError("传入的是 dict，但未找到可用的 DataFrame")
        if hasattr(obj, 'columns'):
            return obj
        raise TypeError("输入既不是 DataFrame 也不是包含 DataFrame 的 dict")

    def get_data_triplet(obj, name="Unknown"):
        df = _resolve_df(obj)
        if 'w' not in df.columns:
            raise KeyError(f"[{name}] 缺少 'w' 列")
        w_vals = np.asarray(df['w'], dtype=float)

        if 's' in df.columns and 'x' in df.columns:
            s_vals = np.asarray(df['s'], dtype=float)
            x_vals = np.asarray(df['x'], dtype=float)
        elif 'x' in df.columns:
            x_vals = np.asarray(df['x'], dtype=float)
            s_vals = x_vals.copy() 
        else:
            raise KeyError(f"[{name}] 缺少坐标列 (s 或 x)")
        
        mask = np.isfinite(s_vals) & np.isfinite(x_vals) & np.isfinite(w_vals)
        s_out, x_out, w_out = s_vals[mask], x_vals[mask], w_vals[mask]
        idx = np.argsort(s_out)
        return s_out[idx], x_out[idx], w_out[idx]

    sF, xF, wF = get_data_triplet(data_FEM, "FEM")
    sR, xR, wR = get_data_triplet(data_RK, "RK")
    sA, xA, wA = get_data_triplet(data_Analysis, "Analysis")

    # --- 2. 插值对齐 (保持不变，用于计算误差) ---
    s_min = max(sF.min(), sR.min(), sA.min())
    s_max = min(sF.max(), sR.max(), sA.max())
    mask_ref = (sF >= s_min) & (sF <= s_max)
    s_ref = sF[mask_ref]
    wF_ref = wF[mask_ref]
    x_ref = xF[mask_ref]

    if len(s_ref) < 2:
        raise ValueError("三组数据公共 s 区间太小，无法对比。")

    fR = interp1d(sR, wR, kind='linear', bounds_error=False, fill_value=np.nan)
    fA = interp1d(sA, wA, kind='linear', bounds_error=False, fill_value=np.nan)
    wR_ref = fR(s_ref)
    wA_ref = fA(s_ref)

    valid = np.isfinite(wF_ref) & np.isfinite(wR_ref) & np.isfinite(wA_ref)
    wF_ref, wR_ref, wA_ref = wF_ref[valid], wR_ref[valid], wA_ref[valid]
    x_ref = x_ref[valid]

    # --- 3. 误差指标计算 (保持不变) ---
    def metrics(y_true, y_pred):
        rmse = np.sqrt(mean_squared_error(y_true, y_pred))
        mae = mean_absolute_error(y_true, y_pred)
        r2 = r2_score(y_true, y_pred)
        rel = np.mean(np.abs(y_pred - y_true) / (np.maximum(np.abs(y_true), 1e-12))) * 100.0
        return rmse, mae, r2, rel

    rmse_R, mae_R, r2_R, rel_R = metrics(wF_ref, wR_ref)
    rmse_A, mae_A, r2_A, rel_A = metrics(wF_ref, wA_ref)
    better_by_rmse = "RK" if rmse_R < rmse_A else "Analysis"
    better_by_r2 = "RK" if r2_R > r2_A else "Analysis"

    # --- 4. 绘图部分 (核心修改在这里！！！) ---
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 子图1：三条曲线
    ax1 = axes[0, 0]

    # === [新增逻辑] 为了视觉对齐，对 Analysis 的 x 坐标进行缩放 ===
    x_max_fem = np.max(xF)       # FEM 的真实末端位置 (e.g., 0.92)
    x_max_ana = np.max(xA)       # Analysis 的原始末端位置 (e.g., 1.0)
    
    # 只要解析解比 FEM 长，就进行压缩
    if x_max_ana > x_max_fem:
        scale_ratio = x_max_fem / x_max_ana
        xA_plot = xA * scale_ratio
        label_suffix = " (Scaled X)" # 在图例里标注一下，说明是缩放过的
    else:
        xA_plot = xA
        label_suffix = ""
    
    # === 绘图 ===
    ax1.plot(xF, wF, 'b-', lw=2.0, label='FEM', alpha=0.9)     
    ax1.plot(xR[::60], wR[::60], 'rx', markersize=10, lw=1.0, label='RK', alpha=0.9)
    # 注意：这里改成了 xA_plot
    ax1.plot(xA_plot[::70], wA[::70], 'g*', markersize=10, lw=1.0, label=f'Analysis{label_suffix}', alpha=0.9)
    
    ax1.set_xlabel('x (Normalized Spatial Coordinate)') # 修改标签，避免误导
    ax1.set_ylabel('w')
    ax1.set_title('Deflection Curves\n(Analysis X-axis scaled to match deformed length)')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.legend(loc='best')

    # 子图2：与 FEM 的绝对差 (保持不变)
    ax2 = axes[0, 1]
    diff_R = np.abs(wR_ref - wF_ref)
    diff_A = np.abs(wA_ref - wF_ref)
    ax2.plot(x_ref, diff_R, 'b-', lw=2, label='|RK - FEM|')
    ax2.plot(x_ref, diff_A, 'r-', lw=2, label='|Analysis - FEM|')
    ax2.fill_between(x_ref, 0, diff_R, color='blue', alpha=0.15)
    ax2.fill_between(x_ref, 0, diff_A, color='red', alpha=0.15)
    ax2.set_xlabel('x'); ax2.set_ylabel('Absolute Difference')
    ax2.set_title('Absolute Error vs FEM\n(Aligned by Material Coordinate s)') 
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.legend(loc='best')

    # 子图3：散点对比 (保持不变)
    ax3 = axes[1, 0]
    ax3.scatter(wF_ref, wR_ref, c='blue', s=1, alpha=1, edgecolors='none', label='RK vs FEM')
    ax3.scatter(wF_ref, wA_ref, c='red', s=1, alpha=1, edgecolors='none', label='Analysis vs FEM')
    wmin = np.nanmin([wF_ref, wR_ref, wA_ref])
    wmax = np.nanmax([wF_ref, wR_ref, wA_ref])
    ax3.plot([wmin, wmax], [wmin, wmax], 'g--', lw=1.5, alpha=0.8, label='y = x')
    ax3.set_xlabel('FEM w')
    ax3.set_ylabel('Predicted w')
    ax3.set_title('Correlation: RK/Analysis vs FEM')
    ax3.grid(True, alpha=0.3, linestyle='--')
    ax3.legend(loc='best', frameon=True, handles=[
            plt.Line2D([0], [0], color='blue', linestyle='-', lw=2, label='RK vs FEM'),
            plt.Line2D([0], [0], color='red', linestyle='-', lw=2, label='Analysis vs FEM'),
            plt.Line2D([0], [0], color='green', linestyle='--', lw=1.5, label='y = x')
        ])

    # 子图4：全局指标 (保持不变)
    ax4 = axes[1, 1]
    names = ['RK', 'Analysis']
    rmses = [rmse_R, rmse_A]
    maes = [mae_R, mae_A]
    r2s = [r2_R, r2_A]
    rels = [rel_R, rel_A]
    xbar = np.arange(len(names))
    width = 0.35
    b1 = ax4.bar(xbar - width/2, rmses, width=width, color='#4C78A8', label='RMSE')
    b2 = ax4.bar(xbar + width/2, maes,  width=width, color='#F58518', label='MAE')
    ax4_t = ax4.twinx()
    b3 = ax4_t.plot(xbar, r2s, 'ks-', lw=2, label='R2')
    ax4.set_xticks(xbar); ax4.set_xticklabels(names)
    ax4.set_ylabel('Error (RMSE/MAE)')
    ax4_t.set_ylabel('R2')
    ax4.set_title('Global Errors (RK/Analysis vs FEM)')
    ax4.grid(True, axis='y', alpha=0.3, linestyle='--')
    y_max = max(max(rmses), max(maes))
    dy = 0.02 * y_max if y_max > 0 else 1e-12
    for i, (rm, ma, r2v, relv) in enumerate(zip(rmses, maes, r2s, rels)):
        ax4.text(i - width/2, rm + dy, f'RMSE={rm:.2e}', ha='center', va='bottom', fontsize=9) 
        ax4.text(i + width/2, ma + dy, f'MAE={ma:.2e}',  ha='center', va='bottom', fontsize=9)
        ax4_t.annotate(f'R2={r2v:.3f}\nRel={relv:.2f}%', (xbar[i], r2v), textcoords='offset points', xytext=(0, 6), ha='center', va='bottom', fontsize=9)
    h1, l1 = ax4.get_legend_handles_labels()
    h2, l2 = ax4_t.get_legend_handles_labels()
    ax4.legend(h1 + h2, l1 + l2, loc='best')

    plt.tight_layout()
    plt.savefig(save_name, dpi=300, bbox_inches='tight')
    plt.show()

    # 文本报告
    print("\n" + "="*80)
    print("THREE-SET COMPARISON REPORT (vs FEM)")
    print("="*80)
    print("Note: In Plot (a), Analysis X-axis is scaled to align with FEM deformed length.")
    print("      Error metrics are calculated based on Material Coordinate 's' (Exact alignment).")
    print(f"Common data points: {len(x_ref)}")
    print("\n[RK vs FEM]")
    print(f"  RMSE = {rmse_R:.6e} | MAE = {mae_R:.6e} | R2 = {r2_R:.6f} | Mean Rel Diff = {rel_R:.3f}%")
    print("\n[Analysis vs FEM]")
    print(f"  RMSE = {rmse_A:.6e} | MAE = {mae_A:.6e} | R2 = {r2_A:.6f} | Mean Rel Diff = {rel_A:.3f}%")
    print("\nFigure saved as:", save_name)
    print("="*80 + "\n")

if __name__ == "__main__":
    data_FEM = opensees_beam.run_beam_opensees(
        L=L,               
        E=E,             
        I=I,           
        A=A,              
        n_elem=100,          
        bc="cantilever",     
        point_loads=[
            {"x": 0.5, "Fy": 0},  
            {"x": 1, "Mz": 0} 
        ],
        compute_N=True,
        sample_n=2000,        
        geometric_nonlinearity=True,  
        dist_loads=[{'x0': 0.0, 'x1': 1.0, 'qy': -10000}],
        verbose=False
    )
    # 确保 RK 代码里的参数也正确 (I = 1.33e-08)
    data_RK = RK.result_df
    analyzer = analysis_solve.AnalysisFunc() 
    data_Analysis = analyzer.generate_situation_4_data(
        M=0,  
        E=E, 
        I=I, 
        l=L, 
        q=10000,
        a=0.5,
        F=0,
        num_points=2000
    )

    compare_three_datasets(data_Analysis=data_Analysis, data_FEM=data_FEM, data_RK=data_RK, save_name='three.png')
