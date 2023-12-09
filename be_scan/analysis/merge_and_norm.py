"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Aggregate raw counts and perform log2-transform and t0 normalization}
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd

def merge_and_norm(sample_sheet, in_ref, 
                   file_dir='',
                   t0='t0', dir_counts='', 
                   out='agg_log2_t0.csv', 
                   save=True, return_df=True,
                   ):
    """
    For a given set of samples and their raw read count files from count_reads,
    aggregate them into a single dataframe, normalize counts to reads per million,
    perform log2 transform, and then normalize to t0. The aggregated raw reads,
    log2 transformed values, and t0 normalized values can be saved to csv.

    Parameters
    ----------
    sample_sheet : str or path
        REQUIRED COLS: 'condition', 'counts_file'
        a sheet with information on sequence id, 
        in_fastq (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_np (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    file_dir : str or path, defaults to ''
        String or path to the directory where all files are found and saved. 
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
    # import sample_sheet and extract information for conditions and associated files
    df_samples = pd.read_csv(sample_sheet)
    dict_counts = dict(zip(df_samples.condition, df_samples.counts_file))

    if 'sgRNA_seq' not in df_ref.columns.tolist():
        raise Exception('in_ref is missing column: sgRNA_seq')
    if t0 not in dict_counts.keys():
        raise Exception ('sample sheet is missing the t0 sample')
    # rearrange the dict samples to place t0 first (for order purposes later)
    list_samples = [t0] + [samp for samp in dict_counts.keys() if samp != t0]
    # aggregate read counts from all samples into df_rawreads
    # also perform log2 norm (brian/broad method; log2(rpm + 1 / total reads))
    df_reads = df_ref.copy()
    for sample in list_samples:
        filepath = file_dir + dict_counts[sample]
        df_temp = pd.read_csv(inpath / filepath, names=['sgRNA_seq', sample])
        # aggregating raw reads
        df_reads = pd.merge(df_reads, df_temp, on='sgRNA_seq')
        # log2 normalization
        total_reads = df_reads[sample].sum()
        df_reads[sample+'_log2'] = df_reads[sample].apply(lambda x: np.log2((x * 1e6 / total_reads) + 1))
        # t0 normalization
        df_reads[sample+'_t0'] = df_reads[sample+'_log2'].sub(df_reads[t0+'_log2'])
    # drop the t0 column since it will be 0
    df_reads.drop(columns=t0+'_log2', inplace=True)

    # export files and return dataframes if necessary
    outpath = path / file_dir
    Path.mkdir(outpath, exist_ok=True)
    # determine which files to export
    if save == True:
        df_reads.to_csv(outpath / out, index=False)
    # determine df to return
    print('Merge and normalize completed')
    if return_df: 
        return df_reads
    