"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Aggregate raw counts and perform log2-transform and t0 normalization}
"""

from pathlib import Path
import pandas as pd

def compare_conds(list_comparisons, in_lfc, in_ref, save=True, out_folder='',
                  out_comps='comparisons.csv', return_df=False):
    """
    Perform pairwise comparisons given a list and export the output to a csv.

    Given a list of comparisons (e.g. treatment vs. control), perform pairwise
    comparisons, generate a dataframe, and export to csv. The list of comparisons
    must be in the format (comparison name, condition 1, condition 2).
    The comparison is performed as (condition 1 - condition 2). Note that this
    can be applied to any format of values, not just averaged condition reps.

    Parameters
    ----------
    list_comparisons : list of tuples in format (name, sample 1, sample 2)
        A list of tuples denoting the comparisons to make, with the comparison
        being sample 1 - sample 2 (e.g. treatment - control). The output column
        headers will be labeled by the comparison name in the tuple.
    in_lfc : str or path
        String or path to the csv file containing the values for comparison.
        The column headers must match the sample names in list_comparisons
    in_ref : str or path
        String or path to the reference file. in_ref must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.
    save : bool, default True
        Whether to save the comparisons dataframe as a csv file.
    out_folder : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_comps : str, default 'comparisons.csv'
        Name of the comparisons csv output file.
    return_df : bool, default False
        Whether to return the comparisons dataframe. The default is False.
    """

    # import files, define variables, check for requirements
    path = Path.cwd()
    df_lfc = pd.read_csv(in_lfc)
    df_comps = pd.read_csv(in_ref)
    # perform treatment vs. control comparison
    for name, treatment, control in list_comparisons:
        df_comps[name] = df_lfc[treatment].sub(df_lfc[control])
    # export files and return dataframes if necessary
    if save:
        outpath = path / out_folder
        Path.mkdir(outpath, exist_ok=True)
        df_comps.to_csv(outpath / out_comps, index=False)
    print('Compare conditions completed')
    if return_df:
        return df_comps
    else:
        return
