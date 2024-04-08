"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Run count_reads, merge_and_norm, average_reps, compare_conds}
"""

import pandas as pd

from be_scan.analysis.count_reads import count_reads
from be_scan.analysis.merge_and_norm import merge_and_norm
from be_scan.analysis.average_reps import average_reps
from be_scan.analysis.compare_conds import compare_conds
from be_scan.analysis.calc_controls import calc_controls

def batch_process(sample_sheet, annotated_lib, comparisons, 
                  neg_ctrl_col, neg_ctrl_conditions, 

    KEY_INTERVAL=(10,80), KEY='CGAAACACC', KEY_REV='GTTTGAGA', dont_trim_G=False,
    file_dir='', controls=['t0'], 
    
    out_dir='', out_counts='counts_library.csv', out_lfc='agg_log2_t0.csv', 
    out_conds='avg_conds.csv', out_comps='conditions.csv', out_stats = 'stats.txt',
    save=True, return_df=False,
    ):
    
    """[Summary]
    Given a set of sgRNA sequences, a sample sheet, and a list of comparisons, 
    complete count_reads, merge_and_norm, average_reps, compare_conds runs.

    Parameters
    ----------
    sample_sheet : str or path
        REQUIRED COLS: 'fastq_file', 'counts_file', 'noncounts_file', 'stats_file'
        a sheet with information on sequence id, 
        in_fastq (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_np (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    annotated_lib : str or path
        String or path to the reference file. annotated_lib must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    comparisons : comparisons .csv in format (name, treatment, control)
        A dataframe denoting the comparisons to make, with the comparison
        being treatment - control. The output column
        headers will be labeled by the name in the dataframe.

    file_dir : str or path, defaults to ''
        String or path to the directory where all files are found and saved. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe

    KEY_INTERVAL : tuple, default (10,80)
        Tuple of (KEY_START, KEY_END) that defines the KEY_REGION. Denotes the
        substring within the read to search for the KEY.
    KEY : str, default 'CGAAACACC'
        Sequence that is expected upstream of the spacer sequence. The default
        is the end of the hU6 promoter.
    KEY_REV : str, default 'GTTTTAGA'
        Sequence that is expected downstream of the spacer sequence. The
        default is the start of the sgRNA scaffold sequence.
    dont_trim_G : bool, default False
        Whether to trim the first G from 21-nt sgRNA sequences to make them 20-nt.

    controls : str, default ['t0']
        Name of the control condition samples in sample_sheet. 
    dir_counts : str, default ''
        Name of the subfolder to find the read count csv files. The default is
        the current working directory.
    out_reads : str, default 'agg_reads.csv'
        Name of the aggregated raw reads csv output file.
    out_log2 : str, default 'agg_log2.csv'
        Name of the aggregated log2 normalized values csv output file.
    out_t0 : str, default 'agg_t0_reps.csv'
        Name of the aggregated t0 normalized values csv output file.

    out_conds : str, default 'agg_t0_conds.csv'
        Name of the averaged replicate values csv output file.

    out_comps : str, default 'agg_comps.csv'
        Name of the comparisons csv output file.
    """

    count_reads_params = {
        'sample_sheet':sample_sheet, 'annotated_lib':annotated_lib, 'file_dir':file_dir, 
        'KEY_INTERVAL':KEY_INTERVAL, 'KEY':KEY, 'KEY_REV':KEY_REV, 'dont_trim_G':dont_trim_G, 
        'out_dir':out_dir, 'out_file':out_counts, 'save':save, 'return_df':return_df, 'save_files':save, 
    }
    count_reads(**count_reads_params)

    merge_and_norm_params = {
        'sample_sheet':sample_sheet, 'counts_library':out_dir+out_counts, 'controls':controls, 
        'out_dir':out_dir, 'out_file':out_lfc, 'save':save, 'return_df':return_df,
    }
    merge_and_norm(**merge_and_norm_params)
    
    average_reps_params = {
        'sample_sheet':sample_sheet, 'log2_subt0':out_dir+out_lfc, 
        'out_dir':out_dir, 'out_file':out_conds, 'save':save, 'return_df':return_df,
    }
    average_reps(**average_reps_params)
    
    compare_conds_params = {
        'comparisons':comparisons, 'avg_conds':out_dir+out_conds, 
        'out_dir':out_dir, 'out_file':out_comps, 'save':save, 'return_df':return_df,
    }
    compare_conds(**compare_conds_params)

    comps = pd.read_csv(comparisons)
    stats_comparisons = comps.name.tolist()
    calc_controls_params = {
        'conditions':out_dir+out_comps, 'stats_comparisons':stats_comparisons, 
        'neg_ctrl_col':neg_ctrl_col, 'neg_ctrl_conditions':neg_ctrl_conditions, 
        'out_dir':out_dir, 'out_file':out_stats, 'save':save, 'return_txt':return_df,
    }
    calc_controls(**calc_controls_params)
