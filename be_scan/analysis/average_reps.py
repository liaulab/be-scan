"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: For a set of conditions, average the replicates and export it to csv}
"""

from pathlib import Path
import warnings

import pandas as pd
import re

def average_reps(sample_sheet, log2_subt0, 
                 
    out_dir='', out_file='avg_conds.csv', 
    save=True, return_df=True,
    ):

    """[Summary]
    Averages the replicates for each condition (e.g. treatment, control) and
    exports the csv file with averaged replicate values for each condition.
    Note that this function is generally used for averaging t0 normalized
    replicate values, but can be used on any set of replicate values
    (e.g. log2 fc values or even raw reads).

    Parameters
    ----------
    sample_sheet : str or path
        REQUIRED COLS: 'condition', 'agg_conditions'
        a sheet with information on sequence id, 
        in_fastq (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_np (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    log2_subt0 : str or path
        String or path to the individual replicate values csv file. The column
        headers for the replicate values must match the keys in dict_conds

    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'avg_conds.csv'
        Name of output dataframe with guides and counts. 
    save : bool, default True
        Whether to save the averaged replicate values as a csv file.
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    """

    # import files, define variables, check for requirements, etc.
    path = Path.cwd()
    df_lfc = pd.read_csv(log2_subt0)
    # import sample_sheet and derive relevant information for conditions and aggregation
    df_samples = pd.read_csv(sample_sheet)
    dict_conds = dict(zip(df_samples.condition, df_samples.agg_conditions))

    # make df to map replicates to condition
    df_map = pd.DataFrame(data=dict_conds.items(), columns=['rep','condition'])

    # generate df to hold the averaged replicates per condition
    for cond in df_map['condition'].unique().tolist():
        # for each unique condition, find the reps
        reps = [x+'_subt0' for x in df_map.loc[df_map['condition'] == cond]['rep'].tolist()]
        # skip averaging for single replicates (otherwise breaks script)
        if len(reps) > 1:
            df_lfc[cond+'_subctrl_avg'] = df_lfc[reps].mean(axis=1)
            df_lfc[cond+'_subctrl_stdev'] = df_lfc[reps].std(axis=1)
        elif len(reps) == 1:
            df_lfc[cond+'_subctrl_avg'] = df_lfc[reps]
            df_lfc[cond+'_subctrl_stdev'] = 0
        else: 
            raise Exception('Error! Replicate number not valid')

    # export files and return dataframes if necessary
    if save:
        outpath = path / out_dir
        Path.mkdir(outpath, exist_ok=True)
        df_lfc.to_csv(outpath / out_file, index=False)
        print('average_reps outputed to', str(outpath / out_file))
    print('Average reps completed')
    if return_df: 
        return df_lfc
