"""
3D Isometric Rendering Script.

Generates 3D isometric stress contours based on OpenSees FEM results.
Features geometric scaling and structured camera viewpoints for academic presentation.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from pathlib import Path
from core import config
from core.opensees_beam import get_fem_result

def render_3d_stress_contour(target_folder="single_runs/3d_renders"):
    # Accepts target_folder to separate single runs from batch runs
    output_dir = Path(f"results/{target_folder}").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    
    p = config.PARAMS.copy()
    case_id = p["CASE_ID"]
    L, b, h, I = p["L"], p["b"], p["h"], p["I"]
    
    print(f"Retrieving FEM physical field data for Case {case_id}...")
    df_fem = get_fem_result(p)
    s, w, theta, M = df_fem['s'].values, df_fem['w'].values, df_fem['theta'].values, df_fem['M'].values
    
    if np.max(np.abs(w)) < 1e-12:
        print("Warning: Deformation is negligible. Consider increasing load parameters in config.py.")
        return

    print("Generating 3D isometric stress rendering...")
    
    Stress_top = (-(M * (h/2)) / I) / 1e6  
    Stress_bot = ((M * (h/2)) / I) / 1e6   
    
    X_top = s - (h/2) * np.sin(theta)
    Z_top = w + (h/2) * np.cos(theta)
    X_bot = s + (h/2) * np.sin(theta)
    Z_bot = w - (h/2) * np.cos(theta)
    
    plt.style.use('default')
    fig = plt.figure(figsize=(14, 8), facecolor='white')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('white')
    ax.axis('off')
    
    vmin, vmax = min(np.min(Stress_top), np.min(Stress_bot)), max(np.max(Stress_top), np.max(Stress_bot))
    # Enforce minimal color gradient for pure moment cases to prevent zero-division errors
    if np.isclose(vmin, vmax):
        vmin, vmax = vmin - 1, vmax + 1
        
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.jet
    
    # 1. Undeformed configuration outline
    X_orig, Y_orig = np.meshgrid(s, np.linspace(-b/2, b/2, 2))
    Z_orig_top = np.full_like(X_orig, h/2)
    ax.plot_surface(X_orig, Y_orig, Z_orig_top, color='#cccccc', alpha=0.2, shade=False)
    ax.plot(s, np.full_like(s, -b/2), np.full_like(s, h/2), color='#888888', lw=1, linestyle='--')
    ax.plot(s, np.full_like(s, b/2), np.full_like(s, h/2), color='#888888', lw=1, linestyle='--')

    # 2. Solid stress contour surfaces
    Y_grid = np.linspace(-b/2, b/2, 5)
    for i in range(len(Y_grid) - 1):
        Y1, Y2 = np.full_like(s, Y_grid[i]), np.full_like(s, Y_grid[i+1])
        ax.plot_surface(np.array([X_top, X_top]), np.array([Y1, Y2]), np.array([Z_top, Z_top]), 
                        facecolors=cmap(norm(np.array([Stress_top, Stress_top]))), shade=False, linewidth=0, antialiased=True)
        ax.plot_surface(np.array([X_bot, X_bot]), np.array([Y1, Y2]), np.array([Z_bot, Z_bot]), 
                        facecolors=cmap(norm(np.array([Stress_bot, Stress_bot]))), shade=False, linewidth=0, antialiased=True)

    # High-contrast boundary edge lines
    ax.plot(X_top, np.full_like(s, -b/2), Z_top, color='black', lw=1.5) 
    ax.plot(X_top, np.full_like(s, b/2), Z_top, color='black', lw=1.5)  
    ax.plot(X_bot, np.full_like(s, -b/2), Z_bot, color='black', lw=1.5) 
    ax.plot(X_bot, np.full_like(s, b/2), Z_bot, color='black', lw=1.5)  

    # ==========================================
    # 3. Camera Position and Clipping Optimization
    # ==========================================
    ax.set_xlim([0, L])
    ax.set_ylim([-b*2, b*2]) 
    
    z_min = min(np.min(Z_bot), -h/2)
    z_max = max(np.max(Z_top), h/2)
    z_range = z_max - z_min
    # Add 20% visual margins to the Z-axis bounding box
    ax.set_zlim([z_min - z_range*0.2, z_max + z_range*0.2])
    
    # Adjust visual box aspect ratio for optimal visualization scale
    visual_z = 10 * (max(z_range, 0.02) / L)
    ax.set_box_aspect((10, 1.5, visual_z)) 
    
    # Isometric camera setup: fixes the fixed-end at the left and free-end at the right
    ax.view_init(elev=20, azim=-75)
    
    # ==========================================
    # 4. Colorbar HUD Overlay
    # ==========================================
    m = cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array([])
    cbar = plt.colorbar(m, ax=ax, shrink=0.45, aspect=20, pad=0.02)
    cbar.set_label('Equivalent Bending Stress (MPa)', color='black', fontsize=12, fontweight='bold')
    
    load_info = ""
    if case_id in [2, 3, 8, 9]: load_info += f"F = {p.get('F', 0)} N "
    elif case_id in[1, 5, 6, 7]: load_info += f"M_e = {p.get('M_e', 0)} N·m "
    elif case_id in [4, 10]: load_info += f"q = {p.get('q', 0)} N/m "
    if case_id in [3, 7, 9]: load_info += f" (at a={p.get('a', 0)} m)"

    clean_load = load_info.replace('\n', '_').replace(' ', '').replace('=', '').replace('·', '').replace('/', '_')
    out_path = output_dir / f"3D_Stress_FEM_Case{case_id}_{clean_load}.jpg"
    
    plt.savefig(out_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Isometric rendering saved to: {out_path}\n")
    plt.close(fig)

if __name__ == "__main__":
    print("-" * 50)
    print("Starting 3D Isometric FEM Render")
    print("-" * 50 + "\n")
    # By default, running this script standalone saves to single_runs
    render_3d_stress_contour(target_folder="single_runs/3d_renders")
