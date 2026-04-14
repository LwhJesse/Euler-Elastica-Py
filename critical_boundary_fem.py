"""
FEM Critical Boundary Analysis Script.

Employs adaptive root-finding algorithms (Brent's method) to locate the precise 
load thresholds where linear Analytical theory deviates from the FEM Corotational 
solution by a specified relative error tolerance (default: 5.0%).
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.optimize import brentq

from core import config
from core.opensees_beam import get_fem_result
from core.analysis_solve import AnalysisFunc

class AdaptiveFemAnalyzer:
    def __init__(self, tolerance_percent=5.0):
        self.x0 = tolerance_percent
        self.p = config.PARAMS.copy()
        self.output_dir = Path("results/boundary_analysis").resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.solver = AnalysisFunc()
        
        # Mapping dict identifying the dominant load variable and if it's a bi-variable case
        self.case_config = {
            1: {'load': 'M_e', 'two_var': False},
            2: {'load': 'F',   'two_var': False},
            3: {'load': 'F',   'two_var': True},
            4: {'load': 'q',   'two_var': False},
            5: {'load': 'M_e', 'two_var': False},
            6: {'load': 'M_e', 'two_var': False},
            7: {'load': 'M_e', 'two_var': True},
            8: {'load': 'F',   'two_var': False},
            9: {'load': 'F',   'two_var': True},
            10: {'load': 'q',  'two_var': False}
        }

    def _get_max_relative_error(self, p_current):
        """Computes the maximum relative deflection error between FEM and Analytical solutions."""
        try:
            df_fem = get_fem_result(p_current)
        except Exception:
            return 999.0 # Treat OpenSees divergence as infinite error

        method_name = f"generate_situation_{p_current['CASE_ID']}_data"
        solver_method = getattr(self.solver, method_name)
        df_ana = solver_method(
            E=p_current["E"], I=p_current["I"], l=p_current["L"], 
            q=p_current.get("q", 0), F=p_current.get("F", 0), 
            M_e=p_current.get("M_e", 0), a=p_current.get("a", 0),
            num_points=1000
        )

        s_fem = df_fem['s'].values if 's' in df_fem.columns else np.linspace(0, p_current["L"], len(df_fem))
        w_fem = df_fem['w'].values
        
        f_ana = interp1d(df_ana['s'].values, df_ana['w'].values, kind='linear', fill_value="extrapolate")
        w_ana_aligned = f_ana(s_fem)

        max_w = np.max(np.abs(w_fem))
        
        # Physical truncation deadband: ignores numerical discretization noise at near-zero loads
        if max_w < 1e-4: 
            return 0.0
            
        if np.isnan(max_w) or max_w > 1.0:
            return 999.0
            
        # Add a 1e-5 (0.01 mm) baseline floor to prevent division by zero in floating point operations
        rel_err = np.max(np.abs(w_fem - w_ana_aligned) / (max_w + 1e-5)) * 100.0
        return rel_err

    def find_critical_load_fast(self, case_id, load_var_name, fixed_a=None, prev_crit=None):
        """Efficient root-finding implementation using Brent's method to identify critical thresholds."""
        p_current = self.p.copy()
        p_current["CASE_ID"] = case_id
        p_current["BC_TYPE"] = "cantilever" if case_id <= 4 else "simply_supported"
        if fixed_a is not None:
            p_current["a"] = fixed_a
            
        def error_residual(load_val):
            # Isolate the target load variable while zeroing out others
            for k in['F', 'M_e', 'q']: p_current[k] = 0.0
            p_current[load_var_name] = load_val
            err = self._get_max_relative_error(p_current)
            return err - self.x0

        # Define load orientation and bracket boundaries
        sign = 1.0 if self.p.get(load_var_name, 1.0) > 0 else -1.0
        L_bound = 0.0 
        
        # Utilize previous critical load for hot-start bounded expansion
        if prev_crit is not None and abs(prev_crit) > 1e-1:
            R_bound = prev_crit * 1.5
        else:
            R_bound = sign * 5000.0
            
        if R_bound == 0: R_bound = sign * 1000.0

        # Dynamically expand right bound until residual becomes positive
        iters = 0
        while error_residual(R_bound) <= 0 and iters < 15:
            R_bound *= 2.0
            iters += 1
            
        # Refine root approximation using Brent's method
        try:
            crit_load = brentq(error_residual, L_bound, R_bound, xtol=1e-3)
            return crit_load
        except Exception as e:
            return R_bound

    def run(self, max_anchor_points=50):
        case_id = self.p["CASE_ID"]
        conf = self.case_config.get(case_id)
        
        if not conf:
            print(f"Error: Unknown CASE_ID = {case_id}")
            return

        load_var_name = conf['load']
        is_two_var = conf['two_var']

        print("-" * 50)
        print(f"FEM Boundary Analysis (CASE_ID = {case_id}, Tolerance = {self.x0}%)")
        print(f"Dominant Var: {load_var_name} | Mode: {'Bi-variable' if is_two_var else 'Single-variable'}")
        print("-" * 50 + "\n")

        if is_two_var:
            L = self.p["L"]
            n_elem = self.p.get("n_elem", 100)
            
            # Align anchor points with exact FEM nodes to eliminate pseudo-errors from load misalignment
            fem_nodes = np.linspace(0.0, L, n_elem + 1)
            valid_nodes = fem_nodes[(fem_nodes >= 0.05 * L) & (fem_nodes <= 0.95 * L)]
            
            if len(valid_nodes) > max_anchor_points:
                step = max(1, len(valid_nodes) // max_anchor_points)
                a_values = valid_nodes[::step]
            else:
                a_values = valid_nodes

            critical_loads =[]
            
            print(f"Scanning load position 'a' across {len(a_values)} aligned nodes...")
            prev_F = None
            for i, a in enumerate(a_values):
                F_crit = self.find_critical_load_fast(case_id, load_var_name, fixed_a=a, prev_crit=prev_F)
                critical_loads.append(F_crit)
                prev_F = F_crit
                if i % 10 == 0:
                    print(f"  -> Progress: {i}/{len(a_values)}, a={a:.3f}m, Critical Load={F_crit:.1f}")

            # Bi-variable plotting routine
            plt.figure(figsize=(10, 6))
            plt.plot(a_values, critical_loads, 'r-', lw=3, label=f'Critical Boundary ($x_0$={self.x0}%)')
            
            y_extreme = min(critical_loads) * 1.2 if min(critical_loads) < 0 else max(critical_loads) * 1.2
            plt.fill_between(a_values, critical_loads, y_extreme, color='#ffcccc', alpha=0.5, label='Nonlinear Zone (Danger)')
            plt.fill_between(a_values, 0, critical_loads, color='#e6f2e6', alpha=0.5, label='Linear Zone (Safe)')
            
            plt.title(f'Case {case_id}: FEM Critical {load_var_name} vs Load Position $a$')
            plt.xlabel('Load Position $a$ (m)')
            plt.ylabel(f'Critical Load {load_var_name}')
            plt.grid(True, linestyle='--')
            plt.legend()
            
            out_path = self.output_dir / f"pure_critical_boundary_fem_case{case_id}.jpg"
            plt.savefig(out_path, dpi=300)
            print(f"\nBoundary diagram generated: {out_path}")
            plt.close()

        else:
            print("Computing single-variable critical threshold and error evolution curve...")
            
            # 1. Probe for the upper bound of the nonlinear spectrum
            sign = 1.0 if self.p.get(load_var_name, -1.0) >= 0 else -1.0
            max_load = sign * 1.0 
            
            p_test = self.p.copy()
            for key in ['F', 'M_e', 'q']: p_test[key] = 0.0
            
            while True:
                p_test[load_var_name] = max_load
                p_test["CASE_ID"] = case_id
                err = self._get_max_relative_error(p_test)
                
                if err >= self.x0 * 1.2 or err == 999.0:
                    break
                if abs(max_load) > 1e7:
                    break
                max_load *= 1.5
                
            # 2. Perform high-density sampling across the identified spectrum
            load_arr = np.linspace(0, max_load, 100)
            err_arr =[]
            
            for val in load_arr:
                if val == 0:
                    err_arr.append(0.0)
                    continue
                p_test[load_var_name] = val
                err = self._get_max_relative_error(p_test)
                err_arr.append(err if err != 999.0 else np.nan)
                
            # 3. Exact computation of the critical intersection point
            exact_crit = self.find_critical_load_fast(case_id, load_var_name, fixed_a=self.p.get("L", 0.3)/2)
            
            # 4. Single-variable plotting routine
            fig, ax = plt.subplots(figsize=(10, 6))
            
            ax.fill_between(load_arr, 0, self.x0, color='#e6f2e6', alpha=0.5, label='Linear Zone (Safe)')
            ax.fill_between(load_arr, self.x0, max(err_arr + [self.x0*1.5]), color='#ffcccc', alpha=0.5, label='Nonlinear Zone (Danger)')
            
            ax.plot(load_arr, err_arr, 'b-', lw=3, label='FEM vs Analytical Error (%)')
            ax.axhline(y=self.x0, color='red', linestyle='--', lw=2, label=f'Tolerance Threshold ({self.x0}%)')
            
            if exact_crit != 0.0 and abs(exact_crit) <= abs(max_load):
                ax.axvline(x=exact_crit, color='green', linestyle=':', lw=2, label=f'Critical Load = {exact_crit:.2f}')
                ax.plot(exact_crit, self.x0, 'ro', markersize=8, zorder=5)
                
                # Intelligent directional arc annotation based on load sign
                x_offset = 60 if exact_crit < 0 else -60
                y_offset = 45
                
                ax.annotate(f'{exact_crit:.2f}', 
                            xy=(exact_crit, 0), xycoords=('data', 'axes fraction'), 
                            xytext=(x_offset, y_offset), textcoords='offset points', 
                            arrowprops=dict(arrowstyle="->", color='green', lw=2.5, connectionstyle="arc3,rad=0.2"),
                            color='green', fontweight='bold', fontsize=12,
                            ha='center', va='center',
                            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="green", alpha=0.9),
                            zorder=10)
                
            ax.set_title(f'Case {case_id}: Relative Error Evolution vs Applied Load ({load_var_name})')
            ax.set_xlabel(f'Applied Load {load_var_name}')
            ax.set_ylabel('Relative Deflection Error (%)')
            ax.grid(True, linestyle='--')
            
            # Reposition legend outside the plot bounding box
            ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5))
            plt.tight_layout() 
            
            out_path = self.output_dir / f"single_var_error_curve_case{case_id}.jpg"
            fig.savefig(out_path, dpi=300, bbox_inches='tight') 
            
            print(f"Computation complete. Critical {load_var_name} = {exact_crit:.2f}")
            print(f"Evolution curve saved to: {out_path}")
            plt.close(fig)

if __name__ == "__main__":
    analyzer = AdaptiveFemAnalyzer(tolerance_percent=5.0)
    analyzer.run(max_anchor_points=150)
