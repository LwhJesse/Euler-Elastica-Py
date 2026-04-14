"""
Batch 3D Isometric Rendering Script.

Iterates through all predefined cases in `core/batch_cases.py` and sequentially 
generates OpenSees 3D isometric stress contours for each.
"""

import matplotlib.pyplot as plt
from core import config
from core.batch_cases import BATCH_TEST_CASES
from run_3d_render import render_3d_stress_contour 

def run_all_3d_renders():
    print("-" * 50)
    print(f"Starting batch execution for {len(BATCH_TEST_CASES)} 3D stress renderings...")
    print("-" * 50 + "\n")

    for case in BATCH_TEST_CASES:
        cid = case['id']
        print(f"Rendering 3D Case: {cid}")
        
        # Override configuration parameters
        config.PARAMS.update(case["p"])
        
        # Re-derive geometry and boundary attributes to prevent parameter mismatch
        b = config.PARAMS.get("b", 0.016)
        h = config.PARAMS.get("h", 0.001)
        config.PARAMS["I"] = (b * h**3) / 12.0
        config.PARAMS["A"] = b * h
        config.PARAMS["BC_TYPE"] = "cantilever" if config.PARAMS["CASE_ID"] <= 4 else "simply_supported"

        try:
            # Call the core rendering logic and explicitly direct to batch_runs folder
            render_3d_stress_contour(target_folder="batch_runs/3d_renders")
            
            # Free memory to prevent figure accumulation and memory leaks
            plt.close('all') 
            
        except Exception as e:
            print(f"Error occurred during 3D rendering of case {cid}: {e}")

    print("\nBatch 3D rendering completed successfully. Check the 'results/batch_runs/3d_renders/' directory.")

if __name__ == "__main__":
    run_all_3d_renders()
