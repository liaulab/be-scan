"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: }
"""

from pathlib import Path
import warnings
import json

import pandas as pd
import re

# For a set of conditions, average the replicates and export it to csv.
def average_reps(sample_sheet, in_lfc, 
                 save=True, out_folder='', out_conds='agg_t0_conds.csv', 
                 return_df=False):
    """
    Averages the replicates for each condition (e.g. treatment, control) and
    exports the csv file with averaged replicate values for each condition.
    Note that this function is generally used for averaging t0 normalized
    replicate values, but can be used on any set of replicate values
    (e.g. log2 fc values or even raw reads).

    Parameters
    ----------
    sample_sheet : 
    sample_sheet : a string of dict in format "{'replicate': 'condition'}"
        Dictionary to map replicates (key) to conditions (value). Replicates
        must match the column headers in the in_lfc file (i.e. sample names),
        otherwise they will be discarded before performing the averaging.
        Example: "{'KN-1': 'unsorted', 'KN-2': 'unsorted', 'KN-3': 'sorted'}"
    in_lfc : str or path
        String or path to the individual replicate values csv file. The column
        headers for the replicate values must match the keys in dict_conds
    save : bool, default True
        Whether to save the averaged replicate values as a csv file.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_conds : str, default 'agg_t0_conds.csv'
        Name of the averaged replicate values csv output file.
    return_df : bool, default False
        Whether to return the averaged reps dataframe The default is False.
    """

    # import files, define variables, check for requirements, etc.
    path = Path.cwd()
    df_lfc = pd.read_csv(in_lfc)
    df_samples = pd.read_csv(sample_sheet)
    dict_conds = dict(zip(df_samples.condition, df_samples.agg_conditions))

    # make df to map replicates to condition
    df_map = pd.DataFrame(data=dict_conds.items(), columns=['rep','condition'])
    # check to make sure replicates are in the input file
    list_reps = df_map['rep'].tolist()
    if not all(rep in list_reps for rep in df_lfc.columns.tolist()):
        list_miss = [rep for rep in list_reps if rep not in df_lfc.columns.tolist()]
        # remove missing replicates from df_map and raise Warning
        df_map = df_map.loc[~df_map['rep'].isin(list_miss)]
        warnings.warn('in_lfc is missing replicates (removed from analysis): ' + str(list_miss))

    # generate df to hold the averaged replicates per condition
    for cond in df_map['condition'].unique().tolist():
        # for each unique condition, find the reps
        reps = df_map.loc[df_map['condition'] == cond]['rep'].tolist()
        # skip averaging for single replicates (otherwise breaks script)
        if len(reps) > 1:
            df_lfc[cond] = df_lfc[reps].mean(axis=1)
        elif len(reps) == 1:
            df_lfc[cond] = df_lfc[reps]
        else:
            raise Exception('Error! Replicate number not valid')

    # export files and return dataframes if necessary
    if save:
        outpath = path / out_folder
        Path.mkdir(outpath, exist_ok=True)
        df_lfc.to_csv(outpath / out_conds, index=False)
    print('Average reps completed')

    if return_df: 
        return df_lfc
