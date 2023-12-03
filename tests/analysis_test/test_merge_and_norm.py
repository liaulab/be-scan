"""
Author: Calvin XiaoYang Hu
Date: 231202

{Description: suit of unit tests for testing functions in be_scan/be_scan/analysis/merge_and_norm}
"""

from be_scan.analysis import merge_and_norm
import os
import pandas as pd

file_dir = "tests/test_data/analysis_data/"

def test_merge_and_norm(): 
    merge_and_norm(sample_sheet   = file_dir + "sample_sheet_batch_count_CBE.csv", 
                   in_ref         = file_dir + "CRAF_and_cntrls_ref_lib.csv",
                   file_dir       = file_dir,
                    )
    df_log2 = pd.read_csv(file_dir + "agg_log2.csv")
    df_reads = pd.read_csv(file_dir + "agg_reads.csv")
    df_t0reps = pd.read_csv(file_dir + "agg_t0_reps.csv")
    assert 't0' in df_log2.columns and 't0' in df_reads.columns
    # this is only t0 so there is no new column added to agg_t0_reps

    # clean up
    os.remove(file_dir + "agg_log2.csv")
    os.remove(file_dir + "agg_reads.csv")
    os.remove(file_dir + "agg_t0_reps.csv")
