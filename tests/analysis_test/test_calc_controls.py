"""
Author: Calvin XiaoYang Hu
Date: 240209

{Description: suit of unit tests for testing functions in be_scan/be_scan/analysis/calc_controls}
"""

from be_scan.analysis import calc_controls
import os
import pandas as pd

file_dir = "tests/test_data/analysis/"

def test_calc_controls(): 
    calc_controls(conditions          = f"{file_dir}compare_conds_sample_out.csv",
                  stats_comparisons   = ["cond1"], 
                  neg_ctrl_col        = "gene", 
                  neg_ctrl_conditions = ["control"],
                  out_dir             = file_dir,
                  )
    with open(f"{file_dir}stats.txt", 'r') as file:
        data = file.read().rstrip()
    assert "For comparison" in data

    # clean up
    os.remove(f"{file_dir}stats.txt")
