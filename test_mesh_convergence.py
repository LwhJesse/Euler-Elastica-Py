"""
FEM Mesh Convergence Verification Script.

Executes a grid convergence study using a "Worst-Case Scenario" approach.
It evaluates a cantilever beam under an extreme concentrated tip load (Case 2) 
that induces severe geometric non-linearity (curvature and large rigid rotations).
Demonstrating convergence under this extreme state guarantees mesh independence 
for all other milder load cases.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from core import config
from core.opensees_beam import get_fem_result

def run_mesh_convergence():
    # Array of element mesh sizes to evaluate
    mesh_sizes =[2, 5, 10, 20, 50, 100, 200, 500]
    max_deflections =[]
    
    # Enforce a worst-case highly nonlinear large-deformation scenario
    p = config.PARAMS.copy()
    p["CASE_ID"] = 2  
    p["F"] = -5000    
    p["BC_TYPE"] = "cantilever"
    
    print("-" * 50)
    print("Initiating FEM mesh convergence study (Worst-case Scenario)...")
    print("-" * 50)
    
    for n in mesh_sizes:
        p["n_elem"] = n
        df_fem = get_fem_result(p)
        max_w = np.max(np.abs(df_fem['w'].values))
        max_deflections.append(max_w)
        print(f"Elements: {n:3d} -> Max Deflection: {max_w:.6f} m")
        
    # Calculate relative error utilizing the finest mesh (N=500) as the baseline truth
    baseline = max_deflections[-1]
    errors =[abs(w - baseline)/baseline * 100 for w in max_deflections]

    # Generate academic convergence plot
    fig, ax1 = plt.subplots(figsize=(8, 5))

    color = 'tab:blue'
    ax1.set_xlabel('Number of Elements (FEM Mesh)')
    ax1.set_ylabel('Max Deflection (m)', color=color)
    ax1.plot(mesh_sizes, max_deflections, 'o-', color=color, linewidth=2, markersize=8)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.grid(True, linestyle='--')

    ax2 = ax1.twinx()  
    color = 'tab:red'
    ax2.set_ylabel('Relative Error to N=500 (%)', color=color)  
    ax2.plot(mesh_sizes, errors, 's--', color=color, linewidth=2, markersize=8)
    ax2.tick_params(axis='y', labelcolor=color)
    
    # Indicate the 1% error stability threshold
    ax2.axhline(1.0, color='red', linestyle=':', alpha=0.5, label='1% Error Threshold')
    
    plt.title('FEM Mesh Convergence Study (Large Deformation)')
    fig.tight_layout()  
    
    # Save to the structured directory
    output_dir = Path("results/mesh_convergence").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / "mesh_convergence.jpg"
    
    plt.savefig(out_path, dpi=300)
    print(f"\nConvergence study complete. Plot saved to: {out_path}")

if __name__ == "__main__":
    run_mesh_convergence()
