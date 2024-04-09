"""
Author: Calvin XiaoYang Hu
Date: 240209

{Description: suit of unit tests for testing functions in be_scan/be_scan/analysis/batch_process}
"""

from be_scan.analysis import batch_process
import os
import pandas as pd

file_dir = "tests/test_data/analysis/"

# def test_batch_process(): 
#     batch_process(sample_sheet   = f"{file_dir}sample_sheet.csv", 
#                   annotated_lib  = f"{file_dir}count_reads_sample_out.csv",
#                   out_dir        = file_dir, 
#                   controls       = ['counts1'],
#                  )
#     df_reads = pd.read_csv(f"{file_dir}agg_log2_t0.csv")
#     assert 'counts1_subctrl_avg' in df_reads.columns

#     # clean up
#     os.remove(f"{file_dir}agg_log2_t0.csv")
