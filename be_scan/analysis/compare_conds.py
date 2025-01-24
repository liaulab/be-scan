"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Compare a conditions across logfc data given a set of comparisons}
"""

from pathlib import Path
import pandas as pd

def compare_conds(comparisons, avg_conds, 
                  
    controls=['t0'], 
    in_dir='', out_dir='', out_file='conditions.csv', 
    save=True, return_df=True, 
    ):
    
    """[Summary]
    Perform pairwise comparisons given a list and export the output to a csv.

    Given a list of comparisons (e.g. treatment vs. control), perform pairwise
    comparisons, generate a dataframe, and export to csv. The list of comparisons
    must be in the format (comparison name, condition 1, condition 2).
    The comparison is performed as (condition 1 - condition 2). Note that this
    can be applied to any format of values, not just averaged condition reps.

    Parameters
    ------------
    comparisons : comparisons .csv in format (name, treatment, control)
        A dataframe denoting the comparisons to make, with the comparison
        being treatment - control. The output column
        headers will be labeled by the name in the dataframe.
    avg_conds : str or path
        String or path to the csv file containing the values for comparison.
        The column headers must match the sample names in comparisons

    in_dir : str or path, defaults to ''
        String or path to the directory where library_counts and sample_sheet is found. 
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'conditions.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe

    Returns
    ------------
    """

    # import files, define variables, check for requirements
    in_path = Path(in_dir)
    out_path = Path(out_dir)
    df_conds = pd.read_csv(in_path / avg_conds)

    if len(controls) > 0: suffix = '_LFCminusControl'
    else: suffix = '_LFC'

    comparisons_df = pd.read_csv(in_path / comparisons)
    comparisons_list = list(comparisons_df.itertuples(index=False, name=None))
    # perform treatment vs. control comparison
    for name, treatment, control in comparisons_list:
        df_conds[name] = df_conds[treatment+suffix+'_avg'].sub(df_conds[control+suffix+'_avg'])

    # export files and return dataframes if necessary
    if save: 
        Path.mkdir(out_path, exist_ok=True)
        df_conds.to_csv(out_path / out_file, index=False)
        print('compare_conds outputed to', str(out_path / out_file))
    print('Compare conditions completed')
    if return_df:
        return df_conds
