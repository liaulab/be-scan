"""
Author: Calvin XiaoYang Hu
Date: 231202

{Description: suit of unit tests for testing functions in be_scan/be_scan/analysis/merge_and_norm}
"""

from be_scan.analysis import compare_conds
import os
import pandas as pd

file_dir = "tests/test_data/analysis_data/"

def test_compare_conds(): 
    compare_conds(in_comparisons = file_dir + "comparisons.csv", 
                  in_conds       = file_dir + "agg_t0_conds_.csv",
                  file_dir       = file_dir,
                  )
    df_comps = pd.read_csv(file_dir + "agg_comps.csv")
    assert 'sorted-unsorted' in df_comps.columns
    # this is only t0 so there is no new column added to agg_t0_reps

    # clean up
    os.remove(file_dir + "agg_comps.csv")
