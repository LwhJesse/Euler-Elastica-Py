"""
Global Configuration Module.

ROLE IN PROJECT:
This file serves as the foundational parameter source for the entire project, 
but its parameters are consumed differently depending on the execution mode:

1. Single-Run Scripts (`run_single.py`, `run_3d_render.py`):
   - Consumes ALL parameters here, including the exact load values (q, F, M_e) 
     and the specified CASE_ID.

2. Batch Scripts (`run_batch.py`, `run_batch_3d.py`):
   - Reads ONLY the base geometry and mesh settings (L, E, b, h, n_elem).
   - The CASE_ID and load values are OVERRIDDEN by `batch_cases.py`.

3. Boundary Search Scripts (`critical_boundary_fem.py`, `run_multiprocess.py`):
   - Reads the base geometry and mesh settings.
   - IGNORES the static load values (q, F, M_e) as the root-finding algorithm 
     dynamically calculates the exact loads that trigger a 5% nonlinear error.
   - Multiprocess mode additionally overrides the CASE_ID to scan all 10 cases.
"""

import matplotlib as mpl
import matplotlib.style as mpl_style

def setup_matplotlib_style():
    """Applies academic formatting to matplotlib globally (Times New Roman & STIX)."""
    _orig_use = mpl_style.use
    def _safe_use(*args, **kwargs):
        _orig_use(*args, **kwargs)
        mpl.rcParams['font.family'] = 'Times New Roman'
        mpl.rcParams['mathtext.fontset'] = 'stix'
    mpl_style.use = _safe_use
    mpl.rcParams['font.family'] = 'Times New Roman'
    mpl.rcParams['mathtext.fontset'] = 'stix'

setup_matplotlib_style()

# =============================================================================
# Global Parameter Dictionary
# =============================================================================
PARAMS = {
    # 1. Active Case Identification (1-10)
    # Applied strictly in single runs and standalone boundary searches.
    "CASE_ID": 6,       # Cases 1-4: Cantilever, Cases 5-10: Simply Supported

    # 2. Base Geometry and Material Properties (Used globally by all scripts)
    "L": 0.3,           # Beam length (m)
    "E": 2e11,          # Young's Modulus (Pa), default: 200 GPa (Steel)
    "b": 0.016,         # Cross-section width (m)
    "h": 0.001,         # Cross-section height (m)

    # 3. Static Load Configurations (Active ONLY in single-run scripts)
    # Note: Downward load directions are considered negative.
    "q": -10000.0,      # Distributed load (N/m) - active for Cases 4, 10
    "F": -3.0,          # Concentrated force (N) - active for Cases 2, 3, 8, 9
    "M_e": -1.0,        # Bending moment (N·m) - active for Cases 1, 5, 6, 7
    "a": 0.015,         # Load action point distance from left support (m)

    # 4. Solver and Discretization Settings
    "n_elem": 200,      # Number of elements for OpenSees FEM
    "sample_n": 1000    # Number of sampling points for analytical functions
}

# ---------------------------------------------------------
# Auto-derived Parameters (Do not modify manually)
# ---------------------------------------------------------
PARAMS["I"] = (PARAMS["b"] * PARAMS["h"]**3) / 12.0  # Moment of inertia (m^4)
PARAMS["A"] = PARAMS["b"] * PARAMS["h"]              # Cross-sectional area (m^2)

# Automatically infer Boundary Condition Type
PARAMS["BC_TYPE"] = "cantilever" if PARAMS["CASE_ID"] <= 4 else "simply_supported"
