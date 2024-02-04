"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Compare a conditions across logfc data given a set of comparisons}
"""

from pathlib import Path
import pandas as pd

def compare_conds(comparisons, 
                  annotated_lib, 
                  out_dir='', out_file='conditions.csv', 
                  save=True, return_df=True, 
                  ):
    """
    Perform pairwise comparisons given a list and export the output to a csv.

    Given a list of comparisons (e.g. treatment vs. control), perform pairwise
    comparisons, generate a dataframe, and export to csv. The list of comparisons
    must be in the format (comparison name, condition 1, condition 2).
    The comparison is performed as (condition 1 - condition 2). Note that this
    can be applied to any format of values, not just averaged condition reps.

    Parameters
    ----------
    comparisons : comparisons .csv in format (name, treatment, control)
        A dataframe denoting the comparisons to make, with the comparison
        being treatment - control. The output column
        headers will be labeled by the name in the dataframe.
    annotated_lib : str or path
        String or path to the csv file containing the values for comparison.
        The column headers must match the sample names in comparisons

    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'conditions.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe
    """

    # import files, define variables, check for requirements
    path = Path.cwd()
    df_conds = pd.read_csv(annotated_lib)

    comparisons_df = pd.read_csv(comparisons)
    comparisons_list = list(comparisons_df.itertuples(index=False, name=None))
    # perform treatment vs. control comparison
    for name, treatment, control in comparisons_list:
        df_conds[name] = df_conds[treatment].sub(df_conds[control])

    # export files and return dataframes if necessary
    if save: 
        outpath = path / out_dir
        Path.mkdir(outpath, exist_ok=True)
        df_conds.to_csv(outpath / out_file, index=False)
        print('compare_conds outputed to', str(outpath / out_file))
    print('Compare conditions completed')
    if return_df:
        return df_conds
