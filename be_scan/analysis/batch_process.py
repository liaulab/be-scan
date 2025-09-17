"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Run count_reads, merge_and_norm, average_reps, compare_conds}
"""

import pandas as pd

from be_scan.analysis.count_reads import count_reads
from be_scan.analysis.log_transform import log_transform
from be_scan.analysis.compare_conds import compare_conds
from be_scan.analysis.calc_controls import calc_controls

def batch_process(
    sample_sheet, 
    annotated_lib, 

    comparisons='', 
    neg_ctrl_col='', neg_ctrl_conditions=[], 
    file_dir='', controls=['t0'], 
    KEY_INTERVAL=(10,80), KEY='CGAAACACC', KEY_REV='GTTTGAGA', dont_trim_G=False,
    sgRNA_seq_col = 'sgRNA_seq', 
    lower_cutoff=0, lower_cutoff_cols=[], 
    out_counts='counts_library.csv', out_lfc='library_LFC.csv', out_comps='conditions.csv', out_stats = 'stats.txt',
    
    stats_comparisons=[], 
    plot_out_type='pdf', save=True, return_df=False,
    ):
    
    """[Summary]
    Given a set of sgRNA sequences, a sample sheet, and a list of comparisons, 
    complete count_reads, merge_and_norm, average_reps, compare_conds runs.

    Parameters
    ------------
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
    neg_ctrl_col : str, defaults to ''
        column of conditions which correspond to normalization control
    neg_ctrl_conditions : list of str, defaults to []
        names of categories of neg_ctrl_col to normalize dataframe

    file_dir : str or path, defaults to ''
        String or path to the directory where all files are found and saved. 
    controls : str, default ['t0']
        Name of the control condition samples in sample_sheet. 

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
    lower_cutoff : float, default 0
        The cutoff for values which do not go into the final dataframe
    lower_cutoff_cols : list of str, default []
        The column names on which to filter for the lower_cutoff

    stats_comparisons : list of str
        a list of columns of the conditions for which to calculate negative controls
    out_dir : str, default ''
        Name of the subfolder to find the read count csv files. The default is
        the current working directory.
    out_counts : str, default 'counts_library.csv'
        Name of the aggregated raw reads csv output file.
    out_lfc : str, default 'library_LFC.csv'
        Name of the aggregated log2 normalized values csv output file.
    out_comps : str, default 'conditions.csv'
        Name of the comparisons csv output file.
    out_stats : str or path, defaults to 'stats.txt'
        Name of output dataframe with guides and counts. 

    plot_out_type : str, optional, defaults to 'pdf'
        file type of figure output for count_reads
    return_df : bool, default False
        Whether or not to return the resulting dataframe and statistics
    save : bool, default True
        Whether or not to save the resulting dataframe

    Returns
    ------------
    """

    count_reads_params = {
        'sample_sheet':sample_sheet, 'annotated_lib':annotated_lib, 'in_dir':file_dir, 
        'KEY_INTERVAL':KEY_INTERVAL, 'KEY':KEY, 'KEY_REV':KEY_REV, 'dont_trim_G':dont_trim_G, 'sgRNA_seq_col':sgRNA_seq_col, 
        'lower_cutoff':lower_cutoff, 'lower_cutoff_cols':lower_cutoff_cols, 
        'out_dir':file_dir, 'out_file':out_counts, 'return_df':return_df, 'plot_out_type':plot_out_type, 'save_files':save, }
    count_reads(**count_reads_params)

    log_transform_params = {
        'sample_sheet':sample_sheet, 'library_counts':out_counts, 'controls':controls, 
        'in_dir':file_dir, 'out_dir':file_dir, 'out_file':out_lfc, 'save':save, 'return_df':return_df, 
    }
    result = log_transform(**log_transform_params)
    
    if len(comparisons) > 1: 
        compare_conds_params = {
            'comparisons':comparisons, 'avg_conds':out_lfc, 'out_file':out_comps, 
            'in_dir':file_dir, 'out_dir':file_dir, 'save':save, 'return_df':return_df, 
            'controls':controls, }
        result = compare_conds(**compare_conds_params)

    if len(neg_ctrl_col) > 0 and len(neg_ctrl_conditions) > 0 and len(stats_comparisons) > 0: 
        calc_controls_params = {
            'conditions':out_comps, 'stats_comparisons':stats_comparisons, 
            'neg_ctrl_col':neg_ctrl_col, 'neg_ctrl_conditions':neg_ctrl_conditions, 
            'in_dir':file_dir, 'out_dir':file_dir, 'out_file':out_stats, 'save':save, 'return_txt':return_df, 
            'controls':controls, }
        calc_controls(**calc_controls_params)

    if return_df: 
        return result
