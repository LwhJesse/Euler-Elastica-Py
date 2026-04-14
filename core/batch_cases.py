"""
Batch Testing Configuration Suite.

ROLE IN PROJECT:
This module stores isolated, predefined parameter sets exclusively used by 
the batch processing scripts (`run_batch.py` and `run_batch_3d.py`). 

These specific load combinations are calibrated to trigger large geometric 
deformations exceeding the 5% linear error threshold, serving as a benchmarking 
suite for the 10 standard beam cases.
"""

BATCH_TEST_CASES =[
    # --- Cantilever Beam Cases ---
    
    # Case 1: End Moment
    {"id": "case1-1", "p": {"CASE_ID": 1, "M_e": -1.0, "F": 0, "q": 0, "a": 0.5}},
    
    # Case 2: End Concentrated Force
    {"id": "case2-1", "p": {"CASE_ID": 2, "F": -3.0, "M_e": 0, "q": 0, "a": 0.5}},
    
    # Case 3: Intermediate Concentrated Force
    {"id": "case3-1", "p": {"CASE_ID": 3, "F": -25.0, "M_e": 0, "q": 0, "a": 0.2}},
    
    # Case 4: Uniformly Distributed Load
    {"id": "case4-1", "p": {"CASE_ID": 4, "q": -25.0, "F": 0, "M_e": 0, "a": 0.5}},
    
    # --- Simply Supported Beam Cases ---
    
    # Case 5: Left End Moment
    {"id": "case5-1", "p": {"CASE_ID": 5, "M_e": -4.0, "F": 0, "q": 0, "a": 0.5}},
    
    # Case 6: Right End Moment (Symmetric to Case 5)
    {"id": "case6-1", "p": {"CASE_ID": 6, "M_e": 4.0, "F": 0, "q": 0, "a": 0.5}},
    
    # Case 7: Intermediate Moment
    {"id": "case7-1", "p": {"CASE_ID": 7, "M_e": -15.0, "a": 0.15, "F": 0, "q": 0}},
    
    # Case 8: Midpoint Concentrated Force
    {"id": "case8-1", "p": {"CASE_ID": 8, "F": -20.0, "M_e": 0, "q": 0, "a": 0.5}},
    
    # Case 9: Intermediate Concentrated Force
    {"id": "case9-1", "p": {"CASE_ID": 9, "F": -25.0, "a": 0.2, "M_e": 0, "q": 0}},
    
    # Case 10: Uniformly Distributed Load
    {"id": "case10-1", "p": {"CASE_ID": 10, "q": -100.0, "F": 0, "M_e": 0, "a": 0.5}},
]
