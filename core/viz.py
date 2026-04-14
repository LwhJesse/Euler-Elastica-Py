"""
Visualization Module.

Handles the generation of 2D comparative plots between FEM, RK, and Analytical solutions.
Enforces an atom-to-atom comparison using Lagrangian coordinates (initial arc length s)
and generates standardized academic error distribution layouts.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
import os
from pathlib import Path

def plot_and_compare(df_fem, df_rk, df_ana, L):
    # Changed output path to structured single_runs directory
    output_dir = Path(__file__).resolve().parent.parent / "results" / "single_runs" / "2d_plots"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    s_base = np.linspace(0, L, 100000)
    
    def align_by_atom(df):
        s_orig = df['s'].values if 's' in df.columns else np.linspace(0, L, len(df))
        aligned = {}
        for col in['theta', 'w', 'x', 'M', 'V', 'N']:
            if col in df.columns:
                f = interp1d(s_orig, df[col].values, kind='linear', fill_value="extrapolate")
                aligned[col] = f(s_base)
        return aligned

    fem = align_by_atom(df_fem)
    rk  = align_by_atom(df_rk)
    ana = align_by_atom(df_ana)

    vars_info =[('w', 'Deflection w (m)'), ('theta', 'Rotation $\\theta$ (rad)'), 
                 ('M', 'Bending Moment M (N·m)'), ('V', 'Shear Force V (N)'), ('N', 'Axial Force N (N)')]

    plt.style.use('default')
    fig1, axes1 = plt.subplots(2, 3, figsize=(18, 10))
    fig1.suptitle('Geometrically Nonlinear Beam Analysis (Lagrangian Description: s-based)', fontsize=16, fontweight='bold')
    
    for i, (var, title) in enumerate(vars_info):
        ax = axes1[i//3, i%3]
        ax.plot(s_base, fem[var], 'k-', lw=2.5, label='FEM (OpenSees)')
        ax.plot(s_base, rk[var],  'r--', lw=2.0, label='RK (Shooting)')
        ax.plot(s_base, ana[var], 'g:',  lw=2.0, label='Analytical')
        ax.set_title(title); ax.set_xlabel('Initial Arc Length s (m)'); ax.grid(True, linestyle='--', alpha=0.3); ax.legend()

    ax_rmse = axes1[1, 2]
    metrics = [v for v, t in vars_info]
    rk_rmse =[np.sqrt(np.mean((rk[v]-fem[v])**2)) / (np.max(fem[v])-np.min(fem[v]) + 1e-9) * 100 for v in metrics]
    ana_rmse = [np.sqrt(np.mean((ana[v]-fem[v])**2)) / (np.max(fem[v])-np.min(fem[v]) + 1e-9) * 100 for v in metrics]
    
    x_pos = np.arange(len(metrics))
    ax_rmse.bar(x_pos - 0.2, rk_rmse, 0.4, label='RK Error %', color='red', alpha=0.6)
    ax_rmse.bar(x_pos + 0.2, ana_rmse, 0.4, label='Ana Error %', color='green', alpha=0.6)
    ax_rmse.set_yscale('log'); ax_rmse.set_xticks(x_pos); ax_rmse.set_xticklabels(metrics)
    ax_rmse.set_ylabel('NRMSE (%) [Log]'); ax_rmse.set_title('Quantized Errors (Lagrangian)'); ax_rmse.legend()

    plt.tight_layout(); fig1.savefig(output_dir / "comparison_result.jpg", dpi=300)

    # Structured layout for absolute error distribution (2 on top, 3 on bottom)
    fig2 = plt.figure(figsize=(18, 10))
    fig2.suptitle('Absolute Error Distribution (Atom-to-Atom)', fontsize=16, fontweight='bold')
    
    import matplotlib.gridspec as gridspec
    gs2 = gridspec.GridSpec(2, 6, figure=fig2)
    
    # Assign specific grid coordinates to achieve centered alignment
    ax2_0 = fig2.add_subplot(gs2[0, 1:3])
    ax2_1 = fig2.add_subplot(gs2[0, 3:5])
    ax2_2 = fig2.add_subplot(gs2[1, 0:2])
    ax2_3 = fig2.add_subplot(gs2[1, 2:4])
    ax2_4 = fig2.add_subplot(gs2[1, 4:6])
    
    axes2_list =[ax2_0, ax2_1, ax2_2, ax2_3, ax2_4]
    
    for i, (var, title) in enumerate(vars_info):
        ax = axes2_list[i]
        d_rk, d_ana = np.abs(rk[var] - fem[var]), np.abs(ana[var] - fem[var])
        
        ax.plot(s_base, d_rk, 'b-', label='|RK - FEM|')
        ax.fill_between(s_base, 0, d_rk, color='blue', alpha=0.15)
        ax.plot(s_base, d_ana, 'r-', label='|Analysis - FEM|')
        ax.fill_between(s_base, 0, d_ana, color='red', alpha=0.15)
        
        ax.set_title(f'Abs Diff of {var}')
        ax.set_xlabel('Initial Arc Length s (m)')
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.legend()

    plt.tight_layout()
    fig2.savefig(output_dir / "error_distribution.jpg", dpi=300)
    
    # Release memory resources to prevent leaks during automated executions
    plt.close(fig1)
    plt.close(fig2)
    print("Lagrangian reconstruction and visualization complete. Results saved.")
