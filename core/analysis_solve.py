"""
Analytical Solution Module for Euler-Bernoulli Beams.

Provides exact analytical solutions (linear small-deflection theory) for 10 standard 
beam configurations. Computes deflection (w), rotation (theta), bending moment (M), 
and shear force (V) along the length of the beam.
"""

import pandas as pd
import numpy as np

class AnalysisFunc:
    def _to_df(self, x, w, theta, m, v, n=None):
        """Helper method to format the analytical results into a standardized DataFrame."""
        if n is None:
            n =[-vi * np.sin(ti) for vi, ti in zip(v, theta)]
        return pd.DataFrame({'s': x, 'theta': theta, 'w': w, 'x': x, 'M': m, 'V': v, 'N': n})

    # =========================================================================
    # Cantilever Beam Cases (1-4)
    # =========================================================================
    
    # Case 1: End Moment (Counter-clockwise M -> Upward deflection)
    def generate_situation_1_data(self, E, I, l, M_e, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        w =[M_e * xi**2 / (2 * E * I) for xi in x]
        theta =[M_e * xi / (E * I) for xi in x]
        return self._to_df(x, w, theta,[M_e]*len(x), [0.0]*len(x))

    # Case 2: End Concentrated Force (Downward F -> Downward deflection)
    def generate_situation_2_data(self, E, I, l, F, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        # Corrected geometric sign convention (3l - xi)
        w =[F * xi**2 / (6 * E * I) * (3 * l - xi) for xi in x]
        theta =[F * xi / (2 * E * I) * (2 * l - xi) for xi in x]
        m_p =[F * (-xi + l) for xi in x]
        v_p = [-F] * len(x)
        return self._to_df(x, w, theta, m_p, v_p)

    # Case 3: Intermediate Concentrated Force (Downward F -> Downward deflection)
    def generate_situation_3_data(self, E, I, l, F, a, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        w, theta, m_p, v_p = [], [], [],[]
        for xi in x:
            if xi <= a:
                w.append(F * xi**2 / (6 * E * I) * (3 * a - xi))
                theta.append(F * xi / (2 * E * I) * (2 * a - xi))
                m_p.append(-F * (xi - a))
                v_p.append(-F)
            else:
                w.append(F * a**2 / (6 * E * I) * (3 * xi - a))
                theta.append(F * a**2 / (2 * E * I))
                m_p.append(0.0)
                v_p.append(0.0)
        return self._to_df(x, w, theta, m_p, v_p)

    # Case 4: Uniformly Distributed Load (Downward q -> Downward deflection)
    def generate_situation_4_data(self, E, I, l, q, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        w =[q * xi**2 / (24 * E * I) * (xi**2 - 4 * l * xi + 6 * l**2) for xi in x]
        theta =[q / (6 * E * I) * (xi**3 - 3 * l * xi**2 + 3 * l**2 * xi) for xi in x]
        m_p = [q / 2.0 * (xi - l)**2 for xi in x]
        v_p = [q * (xi - l) for xi in x]
        return self._to_df(x, w, theta, m_p, v_p)

    # =========================================================================
    # Simply Supported Beam Cases (5-10)
    # =========================================================================
    
    # Case 5: Left End Moment (Counter-clockwise M -> Downward deflection on the left)
    def generate_situation_5_data(self, E, I, l, M_e, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        # Negative sign applied: CCW moment at left support yields downward deflection
        w =[M_e * xi / (6 * E * I * l) * (xi**2 - 3 * l * xi + 2 * l**2) for xi in x]
        theta =[M_e / (6 * E * I * l) * (3 * xi**2 - 6 * l * xi + 2 * l**2) for xi in x]
        m_p =[M_e * (xi / l - 1) for xi in x]
        v_p = [M_e / l] * len(x)
        return self._to_df(x, w, theta, m_p, v_p)

    # Case 6: Right End Moment (Counter-clockwise M -> Upward deflection on the right)
    def generate_situation_6_data(self, E, I, l, M_e, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        # (xi^2 - l^2) yields a negative value, enforcing the upward deflection at the right
        w =[M_e * xi / (6 * E * I * l) * (xi**2 - l**2) for xi in x]
        theta =[M_e / (6 * E * I * l) * (3 * xi**2 - l**2) for xi in x]
        m_p = [M_e * xi / l for xi in x]
        v_p =[M_e / l] * len(x)
        return self._to_df(x, w, theta, m_p, v_p)

    # Case 7: Intermediate Moment
    def generate_situation_7_data(self, E, I, l, M_e, a, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        w, theta, m_p, v_p = [], [], [],[]
        b = l - a
        for xi in x:
            if xi <= a:
                w.append(M_e * xi / (6 * E * I * l) * (xi**2 + 3 * b**2 - l**2))
                theta.append(M_e / (6 * E * I * l) * (3 * xi**2 + 3 * b**2 - l**2))
                m_p.append(M_e * xi / l)
            else:
                w.append(M_e / (6 * E * I * l) * (xi**3 - 3 * l * xi**2 + (3 * a**2 + 2 * l**2) * xi - 3 * a**2 * l))
                theta.append(M_e / (6 * E * I * l) * (3 * xi**2 - 6 * l * xi + 3 * a**2 + 2 * l**2))
                m_p.append(M_e * (xi / l - 1))
            v_p.append(M_e / l)
        return self._to_df(x, w, theta, m_p, v_p)

    # Case 8: Midpoint Concentrated Force (Downward F -> Downward deflection)
    def generate_situation_8_data(self, E, I, l, F, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        w, theta, m_p, v_p = [], [], [],[]
        for xi in x:
            if xi <= l/2:
                w.append(F * xi / (48 * E * I) * (3 * l**2 - 4 * xi**2))
                theta.append(F / (16 * E * I) * (l**2 - 4 * xi**2))
                m_p.append(-F * xi / 2)
                v_p.append(-F / 2)
            else:
                w.append(F * (xi - l) / (48 * E * I) * (4 * (xi - l)**2 - 3 * l**2))
                theta.append(F / (16 * E * I) * (4 * (xi - l)**2 - l**2))
                m_p.append(-F * (l - xi) / 2)
                v_p.append(F / 2)
        return self._to_df(x, w, theta, m_p, v_p)

    # Case 9: Intermediate Concentrated Force (Downward F -> Downward deflection)
    def generate_situation_9_data(self, E, I, l, F, a, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        w, theta, m_p, v_p = [], [], [],[]
        b = l - a
        for xi in x:
            if xi <= a:
                w.append(F * b * xi / (6 * E * I * l) * (l**2 - b**2 - xi**2))
                theta.append(F * b / (6 * E * I * l) * (l**2 - b**2 - 3 * xi**2))
                m_p.append(-F * b * xi / l)
                v_p.append(-F * b / l)
            else:
                w.append(F * a * (l - xi) / (6 * E * I * l) * (l**2 - a**2 - (l-xi)**2))
                theta.append(F * a / (6 * E * I * l) * (3 * (l - xi)**2 + a**2 - l**2))
                m_p.append(-F * a * (l - xi) / l)
                v_p.append(F * a / l)
        return self._to_df(x, w, theta, m_p, v_p)

    # Case 10: Uniformly Distributed Load (Downward q -> Downward deflection)
    def generate_situation_10_data(self, E, I, l, q, **kwargs):
        x = np.linspace(0, l, kwargs.get('num_points', 1000))
        w =[q * xi / (24 * E * I) * (xi**3 - 2 * l * xi**2 + l**3) for xi in x]
        theta =[q / (24 * E * I) * (4 * xi**3 - 6 * l * xi**2 + l**3) for xi in x]
        m_p =[-q * xi / 2.0 * (l - xi) for xi in x]
        v_p =[-q * (l / 2.0 - xi) for xi in x]
        return self._to_df(x, w, theta, m_p, v_p)
