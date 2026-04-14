"""
OpenSeesPy Finite Element Analysis Module.

Constructs and solves a Corotational 2D geometric nonlinear beam model using 
OpenSeesPy. Implements an adaptive stepping algorithm to ensure convergence 
under extreme load conditions.
"""

from typing import Optional, Dict, Any, List, Tuple
import numpy as np
import pandas as pd
import sys
import os

def run_beam_opensees(
    *, L: float, E: float, I: float, A: Optional[float] = None,
    n_elem: int = 50, bc: str = "cantilever",
    point_loads: Optional[List[Dict[str, float]]] = None,
    dist_loads: Optional[List[Dict[str, Any]]] = None,
    compute_N: bool = True, sample_n: Optional[int] = None,
    geometric_nonlinearity: bool = True, max_iterations: int = 50,
    tolerance: float = 1e-8, verbose: bool = False
) -> Dict[str, Any]:
    try:
        import openseespy.opensees as ops
    except ImportError:
        raise ImportError("openseespy not installed.")

    old_stderr = None
    if not verbose:
        old_stderr = sys.stderr
        sys.stderr = open(os.devnull, 'w')

    try:
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)

        n_elem = max(1, int(n_elem))
        node_x = np.linspace(0.0, L, n_elem + 1)
        nnode = node_x.size
        for i, x in enumerate(node_x):
            ops.node(i + 1, float(x), 0.0)

        if bc == "cantilever":
            ops.fix(1, 1, 1, 1)
        elif bc == "simply_supported":
            ops.fix(1, 1, 1, 0)
            ops.fix(nnode, 0, 1, 0)

        mat_tag, sec_tag, transf_tag = 1, 1, 1
        ops.uniaxialMaterial('Elastic', mat_tag, float(E))

        # Apply exact elastic section to eliminate initial moment of inertia errors 
        # that may arise from geometric fiber discretization.
        ops.section('Elastic', sec_tag, float(E), float(A), float(I))

        if geometric_nonlinearity:
            ops.geomTransf('Corotational', transf_tag)
        else:
            ops.geomTransf('Linear', transf_tag)

        n_int = max(3, min(7, int(np.ceil(L / (L / n_elem)))))
        for e in range(n_elem):
            elem_tag = e + 1
            ops.beamIntegration('Lobatto', elem_tag, sec_tag, n_int)
            ops.element('forceBeamColumn', elem_tag, e + 1, e + 2, transf_tag, elem_tag)

        ts_tag, pat_tag = 1, 1
        ops.timeSeries('Linear', ts_tag)
        ops.pattern('Plain', pat_tag, ts_tag)

        if dist_loads:
            for seg in dist_loads:
                qy = float(seg.get('qy', 0.0))
                if abs(qy) > 1e-16:
                    for i in range(1, nnode + 1):
                        trib_len = L / n_elem
                        if i == 1 or i == nnode:
                            trib_len /= 2.0
                        ops.load(i, 0.0, qy * trib_len, 0.0)

        if point_loads:
            for pl in point_loads:
                xp = float(pl.get('x', 0.0))
                node_idx = int(np.argmin(np.abs(node_x - xp)))
                ops.load(node_idx + 1, float(pl.get('Fx', 0.0)), float(pl.get('Fy', 0.0)), float(pl.get('Mz', 0.0)))

        ops.system('BandGeneral')
        ops.numberer('RCM')
        ops.constraints('Plain')

        def solve_adaptive(total_steps_base=20):
            """
            Adaptive solving strategy: Incrementally adjusts the step size and alternates 
            between Newton, ModifiedNewton, and NewtonLineSearch algorithms to maintain 
            convergence in geometrically nonlinear regions.
            """
            attempts = [total_steps_base, 100, 1000]
            
            for n_steps in attempts:
                ops.reset()
                
                ok = True
                step_size = 1.0 / n_steps
                
                ops.integrator('LoadControl', step_size)
                ops.test('NormDispIncr', tolerance, max_iterations)
                ops.algorithm('Newton')
                ops.analysis('Static')
                
                for i in range(n_steps):
                    res = ops.analyze(1)
                    if res != 0:
                        # Fallback to ModifiedNewton if standard Newton fails
                        ops.algorithm('ModifiedNewton')
                        res = ops.analyze(1)
                        if res != 0:
                            # Fallback to Line Search for stiff stabilization
                            ops.algorithm('NewtonLineSearch')
                            res = ops.analyze(1)
                        
                        ops.algorithm('Newton')
                        
                        if res != 0:
                            ok = False
                            break
                
                if ok:
                    return True
            
            return False
        
        # Execute adaptive analysis
        if not solve_adaptive(20 if geometric_nonlinearity else 1):
             raise RuntimeError("OpenSees Analysis failed to converge under adaptive stepping.")

        u_nodes = np.array([ops.nodeDisp(i + 1)[0] for i in range(nnode)])
        w_nodes = np.array([ops.nodeDisp(i + 1)[1] for i in range(nnode)])
        theta_nodes = np.array([ops.nodeDisp(i + 1)[2] for i in range(nnode)])

        x_sections, N_sections, V_sections, M_sections = [], [], [],[]
        for e in range(n_elem):
            forces = ops.eleResponse(e + 1, 'localForces')
            N_sections.append(forces[3]) 
            V_sections.append(-forces[4])  # Ensure shear force direction convention aligns with analytical
            M_sections.append(forces[5])
            x_sections.append(node_x[e + 1])

        if sample_n is None: sample_n = 5 * n_elem
        xs = np.linspace(0.0, L, int(sample_n))

        w_interp = np.interp(xs, node_x, w_nodes)
        theta_interp = np.interp(xs, node_x, theta_nodes)
        Xdef = xs + np.interp(xs, node_x, u_nodes)
        
        final_df = pd.DataFrame({
            "s": xs,
            "theta": theta_interp,
            "M": np.interp(xs, x_sections, M_sections),
            "x": Xdef,
            "w": w_interp,
            "V": np.interp(xs, x_sections, V_sections),
            "N": np.interp(xs, x_sections, N_sections)
        })

        ops.wipe()
        return {"final_df": final_df}

    finally:
        if old_stderr is not None:
            sys.stderr = old_stderr

def get_fem_result(params):
    cid = params.get("CASE_ID", 10)
    L, F, M_e, a, q = params["L"], params.get("F",0), params.get("M_e",0), params.get("a",0.5), params.get("q",0)
    p_loads, d_loads = [],[]
    
    if cid == 1: p_loads = [{'x': L, 'Mz': M_e}]
    elif cid == 2: p_loads = [{'x': L, 'Fy': F}]
    elif cid == 3: p_loads =[{'x': a, 'Fy': F}]
    elif cid == 4: d_loads =[{'x0': 0.0, 'x1': L, 'qy': q}]
    elif cid == 5: p_loads = [{'x': 0.0, 'Mz': M_e}]
    elif cid == 6: p_loads = [{'x': L, 'Mz': M_e}] 
    elif cid == 7: p_loads = [{'x': a, 'Mz': M_e}]
    elif cid == 8: p_loads = [{'x': L/2.0, 'Fy': F}]
    elif cid == 9: p_loads = [{'x': a, 'Fy': F}]
    elif cid == 10: d_loads =[{'x0': 0.0, 'x1': L, 'qy': q}]
    
    res = run_beam_opensees(
        L=L, E=params["E"], I=params["I"], A=params["A"],  
        n_elem=params["n_elem"], bc=params.get("BC_TYPE", "simply_supported"),
        point_loads=p_loads, dist_loads=d_loads,
        geometric_nonlinearity=True, verbose=False
    )
    return res["final_df"][['s', 'theta', 'w', 'x', 'M', 'V', 'N']]
