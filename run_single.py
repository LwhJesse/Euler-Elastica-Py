"""
Single Case Execution Script.

Executes a single beam load case based on the parameters defined in `core/config.py`.
It computes the FEM (OpenSees), RK (Shooting), and Analytical solutions, 
and generates 2D comparison plots.
"""

from core import config
from core.opensees_beam import get_fem_result
from core.RK import get_rk_result
from core.analysis_solve import AnalysisFunc
from core import viz

def main():
    p = config.PARAMS
    cid = p["CASE_ID"]
    
    print("-" * 50)
    print(f"Executing Case ID: {cid} ({p['BC_TYPE']})")
    print("-" * 50)
    
    print("Running OpenSees FEM analysis...")
    df_fem = get_fem_result(p)
    
    print("Running RK Shooting method...")
    df_rk = get_rk_result(p)
    
    print("Computing Analytical solution...")
    solver = AnalysisFunc()
    method_name = f"generate_situation_{cid}_data"
    
    solver_method = getattr(solver, method_name)
    df_ana = solver_method(
        E=p["E"], I=p["I"], l=p["L"], 
        q=p["q"], F=p["F"], M_e=(p["M_e"]), a=p["a"],
        num_points=p["sample_n"]
    )

    print("Computation completed. Generating plots...")
    viz.plot_and_compare(df_fem, df_rk, df_ana, p["L"])

if __name__ == "__main__":
    main()
