"""
Author: Calvin XiaoYang Hu
Date: 231202

{Description: suit of unit tests for testing functions in be_scan/be_scan/analysis/average_reps}
"""

from be_scan.analysis import average_reps
import os
import pandas as pd

file_dir = "tests/test_data/analysis/"

def test_average_reps(): 
    average_reps(sample_sheet   = f"{file_dir}sample_sheet.csv", 
                 log2_subt0     = f"{file_dir}merge_and_norm_sample_out.csv", 
                 out_dir        = file_dir,
                )
    df_conds = pd.read_csv(f"{file_dir}avg_conds.csv")
    assert 'cond1' in df_conds.columns
    assert 'cond1_stdev' in df_conds.columns

    # clean up
    os.remove(f"{file_dir}avg_conds.csv")
