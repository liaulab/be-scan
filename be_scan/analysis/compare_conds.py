"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Compare a conditions across logfc data given a set of comparisons}
"""

from pathlib import Path
import pandas as pd

def compare_conds(in_comparisons, 
                  in_conds, 
                  file_dir='',
                  out_comps='agg_comps.csv', 
                  save=True, return_df=False):
    """
    Perform pairwise comparisons given a list and export the output to a csv.

    Given a list of comparisons (e.g. treatment vs. control), perform pairwise
    comparisons, generate a dataframe, and export to csv. The list of comparisons
    must be in the format (comparison name, condition 1, condition 2).
    The comparison is performed as (condition 1 - condition 2). Note that this
    can be applied to any format of values, not just averaged condition reps.

    Parameters
    ----------
    in_comparisons : in_comparisons .csv in format (name, treatment, control)
        A dataframe denoting the comparisons to make, with the comparison
        being treatment - control. The output column
        headers will be labeled by the name in the dataframe.
    in_conds : str or path
        String or path to the csv file containing the values for comparison.
        The column headers must match the sample names in in_comparisons

    file_dir : str, default ''
        Name of the subfolder to save output files. The default is the current
        working directory.
    out_comps : str, default 'agg_comps.csv'
        Name of the comparisons csv output file.
    return_df : bool, default False
        Whether to return the comparisons dataframe. The default is False.
    """

    # import files, define variables, check for requirements
    path = Path.cwd()
    df_conds = pd.read_csv(in_conds)

    comparisons_df = pd.read_csv(in_comparisons)
    comparisons_list = list(comparisons_df.itertuples(index=False, name=None))
    # perform treatment vs. control comparison
    for name, treatment, control in comparisons_list:
        df_conds[name] = df_conds[treatment].sub(df_conds[control])

    # export files and return dataframes if necessary
    if save: 
        outpath = path / file_dir
        Path.mkdir(outpath, exist_ok=True)
        df_conds.to_csv(outpath / out_comps, index=False)
    print('Compare conditions completed')
    if return_df:
        return df_conds
