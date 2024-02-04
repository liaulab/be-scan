"""
Author: Calvin XiaoYang Hu
Date: 231202

{Description: suit of unit tests for testing functions in be_scan/be_scan/analysis/average_reps}
"""

from be_scan.analysis import average_reps
import os
import pandas as pd

file_dir = "tests/test_data/analysis_data/"

def test_average_reps(): 
    average_reps(sample_sheet   = file_dir + "sample_sheet_batch_count_CBE.csv", 
                 annotated_lib  = file_dir + "CRAF_and_cntrls_ref_lib_agglog2.csv", 
                 out_dir        = file_dir,
                )
    df_conds = pd.read_csv(file_dir + "avg_conds.csv")
    assert 't0_t0' in df_conds.columns

    # clean up
    os.remove(file_dir + "avg_conds.csv")
