"""
Nonlinear Elastica Solver using the Runge-Kutta Shooting Method.

Formulates the exact geometrical nonlinear governing equations of an inextensible 
beam (Euler Elastica). Uses SciPy's solve_ivp (RK45) for forward integration and 
least_squares (Levenberg-Marquardt) to resolve the two-point boundary value problem.
"""

import scipy.integrate as sp
import pandas as pd 
import numpy as np
import math
from scipy.optimize import least_squares

class Event:
    def __init__(self, x_target):
        self.x_target = x_target
        self.terminal = True
        self.direction = 0
    def __call__(self, s, y):
        return s - self.x_target

def integrate_with_u0(u0, load_func, special_point, max_step, start_load, bc_type, L, E, I):
    M_ext_0 = start_load.get('Mz', 0.0)
    
    if bc_type == 'cantilever':
        y0 = [0.0, u0[0], 0.0, 0.0, u0[1], u0[2]]
    else:
        y0 = [u0[0], -M_ext_0, 0.0, 0.0, u0[1], u0[2]]

    def f(s, y):
        theta, M, x, w, V, N = y
        curv = M / (E * I)
        return[curv, V, math.cos(theta), math.sin(theta), 
                N*curv + load_func(s)*math.cos(theta), 
                -V*curv - load_func(s)*math.sin(theta)]

    sorted_sp = sorted(special_point, key=lambda p: p['x'])
    y_sum, t_sum =[],[]
    y_curr = np.array(y0, dtype=float)
    s_start = 0.0
    
    # Dynamic integration step-size matching to prevent missing data points
    # in scenarios with highly localized applied forces (extremely small 'a').
    dynamic_max_step = min(max_step, L / 1000.0)
    
    for i in range(len(sorted_sp) + 1):
        s_end = sorted_sp[i]['x'] if i < len(sorted_sp) else L
        if s_end > s_start:
            sol = sp.solve_ivp(fun=f, t_span=(s_start, s_end), y0=y_curr, 
                               method='RK45', dense_output=True, 
                               max_step=dynamic_max_step, rtol=1e-10, atol=1e-12)
            start_idx = 1 if len(t_sum) > 0 else 0
            for u in range(start_idx, len(sol.t)):
                y_sum.append(sol.y[:, u])
                t_sum.append(sol.t[u])
            y_curr = sol.y[:, -1].copy()
            s_start = s_end
        elif len(t_sum) == 0:
            y_sum.append(y_curr.copy())
            t_sum.append(s_start)
            
        if i < len(sorted_sp):
            f_y = sorted_sp[i].get('force', 0.0)
            y_curr[5] += -f_y * math.sin(y_curr[0])
            y_curr[4] += f_y * math.cos(y_curr[0])
            y_curr[1] -= sorted_sp[i].get('torque', 0.0)
            
    return np.array(t_sum), np.array(y_sum).T

def shooting_solve(start_load, end_load, load_func, special_point, max_step, bc_type, L, E, I):
    Fx = float(end_load.get("Fx", 0.0))
    Fy = float(end_load.get("Fy", 0.0))
    Mz = float(end_load.get("Mz", 0.0))
    M_ext_0 = float(start_load.get("Mz", 0.0))
    
    q_est = load_func(L/2.0)
    
    if bc_type == "cantilever":
        V0_est = Fy - q_est * L - sum(p.get('force', 0.0) for p in special_point)
        M0_est = Mz - Fy * L - q_est * (L**2) / 2.0
        for p in special_point:
            M0_est -= (p.get('torque', 0.0) + p.get('force', 0.0) * p['x'])
        # Strictly enforce moment sign convention for static equilibrium
        M0_est = -M0_est
        u0_guess = np.array([M0_est, V0_est, Fx])
    else:
        moment_balance = Mz + M_ext_0 
        moment_balance += sum(p.get('torque', 0.0) for p in special_point)
        moment_balance -= sum(p.get('force', 0.0) * (L - p['x']) for p in special_point)
        moment_balance -= q_est * (L**2) / 2.0
        
        V0_est = moment_balance / L
        theta0_est = (M_ext_0 * L) / (3.0 * E * I)
        u0_guess = np.array([theta0_est, V0_est, Fx])

    def residual(u0):
        ts, Ys = integrate_with_u0(u0, load_func, special_point, max_step, start_load, bc_type, L, E, I)
        M_e, V_e, N_e, th_e, w_e = Ys[1,-1], Ys[4,-1], Ys[5,-1], Ys[0,-1], Ys[3,-1]
        
        He = N_e * np.cos(th_e) + V_e * np.sin(th_e)
        Ve = N_e * np.sin(th_e) - V_e * np.cos(th_e)
        
        if bc_type == 'cantilever':
            return np.array([(M_e - Mz)/1000, (He - Fx), (Ve - Fy)])
        else:
            return np.array([(M_e - Mz)/1000, w_e * 100, (He - Fx)])

    # Utilize the Levenberg-Marquardt (lm) algorithm for robust convergence
    # in highly nonlinear regions, mitigating Jacobian matrix singularity issues.
    sol = least_squares(residual, u0_guess, method='lm', ftol=1e-12, xtol=1e-12, max_nfev=10000)
    
    if not sol.success and sol.status <= 0:
        import logging
        logging.warning(f"Shooting method convergence issue: {sol.message}")
        
    u0 = sol.x
    ts, Ys = integrate_with_u0(u0, load_func, special_point, max_step, start_load, bc_type, L, E, I)
    
    df = pd.DataFrame({
        "s": ts, "theta": Ys[0,:], "M": Ys[1,:],
        "x": Ys[2,:], "w": Ys[3,:], "V": Ys[4,:], "N": Ys[5,:]
    })
    return df

def get_rk_result(params):
    cid = params.get("CASE_ID", 10)
    L, E, I = params["L"], params["E"], params["I"]
    F_val, M_e_val, a_val = params.get("F", 0.0), params.get("M_e", 0.0), params.get("a", L/2)
    
    q_val = 0.0
    start_load, end_load, sp_points = {}, {}, []
    
    if cid == 1: end_load["Mz"] = M_e_val
    elif cid == 2: end_load["Fy"] = F_val
    elif cid == 3: sp_points.append({"x": a_val, "force": F_val, "torque": 0.0})
    elif cid == 4: q_val = params["q"]
    elif cid == 5: start_load["Mz"] = M_e_val 
    elif cid == 6: end_load["Mz"] = M_e_val
    elif cid == 7: sp_points.append({"x": a_val, "force": 0.0, "torque": M_e_val})
    elif cid == 8: sp_points.append({"x": L/2.0, "force": F_val, "torque": 0.0})
    elif cid == 9: sp_points.append({"x": a_val, "force": F_val, "torque": 0.0})
    elif cid == 10: q_val = params["q"]

    df = shooting_solve(
        start_load, end_load, lambda s: q_val, sp_points, 
        0.0005, params.get("BC_TYPE", "simply_supported"), L, E, I
    )
    return df[['s', 'theta', 'w', 'x', 'M', 'V', 'N']]
