"""
Batch Execution Script.

uterates through the predefined benchmark cases in `core/batch_cases.py`.
For each case, it overrides the config parameters, runs the single-case solver,
and organizes the generated output plots into a dedicated batch directory.
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
from core import config
import run_single
from core.batch_cases import BATCH_TEST_CASES

# 1. Define main output directory
output_dir = Path("results").resolve()
output_dir.mkdir(parents=True, exist_ok=True)

print(f"Starting batch execution for {len(BATCH_TEST_CASES)} predefined cases...\n")

for case in BATCH_TEST_CASES:
    cid = case['id']
    print("-" * 50)
    print(f"Processing Case: {cid}")
    
    # 2. Dynamically override parameters in config.py
    config.PARAMS.update(case["p"])
    
    # Re-derive geometric and boundary properties based on new parameters
    b = config.PARAMS.get("b", 0.02)
    h = config.PARAMS.get("h", 0.02)
    config.PARAMS["I"] = (b * h**3) / 12.0
    config.PARAMS["A"] = b * h
    config.PARAMS["BC_TYPE"] = "cantilever" if config.PARAMS["CASE_ID"] <= 4 else "simply_supported"

    try:
        # 3. Execute main routine (which outputs to 'single_runs/2d_plots' by default)
        run_single.main()  
        
        # 4. Relocate generated plots to a dedicated batch directory
        single_dir = output_dir / "single_runs" / "2d_plots"
        batch_2d_dir = output_dir / "batch_runs" / "2d_plots"
        batch_2d_dir.mkdir(parents=True, exist_ok=True)
        
        orig_comp = single_dir / "comparison_result.jpg"
        orig_err = single_dir / "error_distribution.jpg"
        
        new_comp = batch_2d_dir / f"{cid}_comparison.jpg"
        new_err = batch_2d_dir / f"{cid}_error.jpg"
        
        if new_comp.exists(): new_comp.unlink()
        if new_err.exists(): new_err.unlink()
        
        if orig_comp.exists(): 
            orig_comp.rename(new_comp)
        if orig_err.exists(): 
            orig_err.rename(new_err)
            
        # Free memory to prevent memory leaks during high-res batch rendering
        plt.close('all')
        
    except Exception as e:
        print(f"Error occurred in case {cid}: {e}\n")

print("\nBatch execution completed successfully. Check the 'results/batch_runs/2d_plots/' directory.")
