"""
Author: Calvin XiaoYang Hu
Date: 231202

{Description: suit of unit tests for testing functions in be_scan/be_scan/analysis/merge_and_norm}
"""

from be_scan.analysis import merge_and_norm
import os
import pandas as pd

file_dir = "tests/test_data/analysis/"

def test_merge_and_norm(): 
    merge_and_norm(sample_sheet   = f"{file_dir}sample_sheet.csv", 
                   counts_library = f"{file_dir}count_reads_sample_out.csv",
                   out_dir        = file_dir, 
                   controls       = ['counts1'],
                  )
    df_reads = pd.read_csv(f"{file_dir}agg_log2_t0.csv")
    assert 'counts1_subt0' in df_reads.columns
    assert not ('counts1_log2' in df_reads.columns)

    # clean up
    os.remove(f"{file_dir}agg_log2_t0.csv")
