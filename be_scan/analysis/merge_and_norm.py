"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: }
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import json
import re

# Aggregate raw counts and perform log2-transform and t0 normalization.
def merge_and_norm(sample_sheet, in_ref, 
                   t0='t0', dir_counts='', 
                   file_dir='',
                   save='all',
                   out_folder='', out_reads='agg_reads.csv', out_log2='agg_log2.csv', out_t0='agg_t0_reps.csv', 
                   return_df=None):
    """
    For a given set of samples and their raw read count files from count_reads,
    aggregate them into a single dataframe, normalize counts to reads per million,
    perform log2 transform, and then normalize to t0. The aggregated raw reads,
    log2 transformed values, and t0 normalized values can be saved to csv.

    Parameters
    ----------
    sample_sheet : 
    sample_sheet : a string of dict in format "{'sample name': 'file name'}"
        Dictionary to map sample names (key) to read count file names (value),
        as key,value (e.g. "{'KN-0': 'KN-0_counts.csv'}"). Must include the t0
        sample in the dict. 
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    t0 : str, default 't0'
        Name of the t0 sample in dict_counts. If you have multiple t0 samples
        (e.g. paired t0s for specific samples), then you will need to run this
        function separately for each set of samples with their appropriate t0.
    dir_counts : str, default ''
        Name of the subfolder to find the read count csv files. The default is
        the current working directory.
    save : {'all', None, ['reads', 'log2', 't0']}, default 'all'
        Choose files for export to csv. The default is 'all', which is the
        aggregated read counts ('reads'), log2 normalized values ('log2'), and
        t0 normalized values ('t0'). You may also enter any combination of
        'reads', 'log2', 't0') as a list of strings to choose which ones to
        save ('all' is equivalent to a list of all three). None will not export
        any files to csv.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_reads : str, default 'agg_reads.csv'
        Name of the aggregated raw reads csv output file.
    out_log2 : str, default 'agg_log2.csv'
        Name of the aggregated log2 normalized values csv output file.
    out_t0 : str, default 'agg_t0_reps.csv'
        Name of the aggregated t0 normalized values csv output file.
    return_df : {None, 'reads', 'log2', 't0'}, default None
        Whether to return a dataframe at function end. The default is None,
        which returns nothing. However, you can return the reads, log2 norm,
        or t0 norm values by calling 'reads', 'log2', or't0', respectively.
    """

    # import reference file, define variables, check for requirements
    path = Path.cwd()
    inpath = path / dir_counts
    df_ref = pd.read_csv(in_ref)

    df_samples = pd.read_csv(sample_sheet)
    dict_counts = dict(zip(df_samples.condition, df_samples.counts_file))

    if 'sgRNA_seq' not in df_ref.columns.tolist():
        raise Exception('in_ref is missing column: sgRNA_seq')
    if t0 not in dict_counts.keys():
        raise Exception ('dict_counts is missing the t0 sample')
    # rearrange the dict samples to place t0 first (for order purposes later)
    list_samples = [t0] + [samp for samp in dict_counts.keys() if samp != t0]
    # aggregate read counts from all samples into df_rawreads
    # also perform log2 norm (brian/broad method; log2(rpm + 1 / total reads))
    df_reads, df_log2, df_t0 = df_ref.copy(), df_ref.copy(), df_ref.copy()
    for sample in list_samples:
        filepath = file_dir + dict_counts[sample]
        df_temp = pd.read_csv(inpath / filepath, names=['sgRNA_seq', sample])
        # aggregating raw reads
        df_reads = pd.merge(df_reads, df_temp, on='sgRNA_seq')
        # log2 normalization
        total_reads = df_reads[sample].sum()
        df_log2[sample] = df_reads[sample].apply(lambda x: np.log2((x * 1000000 / total_reads) + 1))
        # t0 normalization
        df_t0[sample] = df_log2[sample].sub(df_log2[t0])
    # drop the t0 column since it will be 0
    df_t0.drop(columns=t0, inplace=True)

    # export files and return dataframes if necessary
    outpath = path / out_folder
    Path.mkdir(outpath, exist_ok=True)
    # dictionary to map kws to dfs and output file names
    dict_df = {'reads': (df_reads, out_reads), 'log2': (df_log2, out_log2), 't0': (df_t0, out_t0)}
    # determine which files to export
    if save == 'all':
        save = ['reads','log2','t0']
    if isinstance(save, list):
        for key in save:
            dict_df[key][0].to_csv(outpath / dict_df[key][1], index=False)
    elif save is None:
        pass
    else:
        warnings.warn('Invalid value for save. No files exported')
    # determine df to return
    print('Merge and normalize completed')
    if return_df in dict_df.keys():
        return dict_df[return_df][0]
    elif return_df is None:
        return
    else:
        print('Invalid value for return_df. No dataframe returned')
        return
    