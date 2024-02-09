"""
Author: Calvin XiaoYang Hu
Date: 231202

{Description: suit of unit tests for testing functions in be_scan/be_scan/analysis/compare_conds}
"""

from be_scan.analysis import compare_conds
import os
import pandas as pd

file_dir = "tests/test_data/analysis/"

def test_compare_conds(): 
    compare_conds(comparisons    = f"{file_dir}comparisons.csv", 
                  annotated_lib  = f"{file_dir}average_reps_sample_out.csv",
                  out_dir        = file_dir, 
                  )
    df_comps = pd.read_csv(f"{file_dir}conditions.csv")
    assert 'cond1-cond1' in df_comps.columns

    # clean up
    os.remove(f"{file_dir}conditions.csv")
