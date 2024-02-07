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

def merge_and_norm(sample_sheet, annotated_lib, 
                   
    controls=['t0'], out_dir='', out_file='agg_log2_t0.csv', 
    save=True, return_df=True,
    ):
    
    """[Summary]
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
    annotated_lib : str or path
        String or path to the reference file. annotated_lib must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.

    controls : str, default ['t0']
        Name of the control condition samples in sample_sheet. 
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'agg_log2_t0.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe
    """

    # import reference file, define variables, check for requirements
    path = Path.cwd()
    df_ref = pd.read_csv(annotated_lib)
    # import sample_sheet and extract information for conditions and associated files
    df_samples = pd.read_csv(sample_sheet)
    dict_counts = dict(zip(df_samples.condition, df_samples.counts_file))

    if 'sgRNA_seq' not in df_ref.columns.tolist():
        raise Exception('annotated_lib is missing column: sgRNA_seq')
    for t0 in controls: 
        if t0 not in dict_counts.keys(): 
            raise Exception ("sample sheet is missing the {0} sample".format(t0))
        # rearrange the dict samples to place t0 first (for order purposes later)
        list_samples = [t0] + [samp for samp in dict_counts.keys() if samp != t0]
        # aggregate read counts from all samples into df_rawreads
        # also perform log2 norm (brian/broad method; log2(rpm + 1 / total reads))
        df_reads = df_ref.copy()
        for sample in list_samples:
            # log2 normalization
            total_reads = pd.to_numeric(df_reads[sample]).sum()
            df_reads[sample+'_log2'] = df_reads[sample].apply(lambda x: np.log2(((float(x) * 1e6 + 1)/ total_reads)))
            # t0 normalization
            df_reads[sample+'_subt0'] = df_reads[sample+'_log2'].sub(df_reads[t0+'_log2'])
        # drop the t0 column since it will be 0
        df_reads.drop(columns=t0+'_log2', inplace=True)

    # export files and return dataframes if necessary
    if save == True: 
        outpath = path / out_dir
        Path.mkdir(outpath, exist_ok=True)
        df_reads.to_csv(outpath / out_file, index=False)
        print('merge_and_norm outputed to', str(outpath / out_file))
    # determine df to return
    print('Merge and normalize completed')
    if return_df: 
        return df_reads
    