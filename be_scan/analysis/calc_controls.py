"""
Author: Calvin XiaoYang Hu
Date: 231128

{Description: Calculate the control stats from the generated dataset}
"""

from pathlib import Path
import pandas as pd
from _controls_ import *

def calc_controls(conditions, stats_comparisons, 
                  neg_ctrl_col, neg_ctrl_conditions,

    in_dir='', out_dir='', out_file='stats.txt', 
    save=True, return_txt=False,
    ):

    """[Summary]
    Calculate and output the normalization controls of the dataframe conditions. 

    Parameters
    ----------
    conditions : str or path
        String or path to the csv file containing the values for comparison.
        The column headers must match the sample names in stats_comparisons
    stats_comparisons : list of str
        a list of columns of the conditions for which to calculate negative controls
    neg_ctrl_col : str
        column of conditions which correspond to normalization control
    neg_ctrl_conditions : list of str
        names of categories of neg_ctrl_col to normalize dataframe
    
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'stats.txt'
        Name of output stats file
    save : bool, default True
        Whether or not to save the resulting dataframe
    return_txt : bool, default False
        Whether or not to print the resulting statistics
    """
    in_path = Path(in_dir)
    out_path = Path(out_dir)
    
    df_conds = pd.read_csv(in_path / conditions)
    stats_comparisons = [x+'_LFCminusControl_avg' for x in stats_comparisons]

    # calculate negative control stats
    _, list_negctrlstats, _ = calc_neg_ctrls(df_conds, stats_comparisons, 
                                             neg_ctrl_col, neg_ctrl_conditions)

    Path.mkdir(out_path, exist_ok=True)
    if save: 
        f = open(out_path / out_file, "w")
        for comp, mean, stdev, top, bottom in list_negctrlstats: 
            f.write("For comparison {0}\n".format(comp))
            f.write("Mean is {0}\n".format(mean))
            f.write("Standard deviation is {0}\n".format(stdev))
            f.write("Mean +- 2 standard deviations {0} {1}\n".format(top, bottom))
        f.close()
        print('calc_controls outputed to', str(out_path / out_file))
    if return_txt: 
        for comp, mean, stdev, top, bottom in list_negctrlstats: 
            print("For comparison {0}\n".format(comp))
            print("Mean is {0}\n".format(mean))
            print("Standard deviation is {0}\n".format(stdev))
            print("Mean +- 2 standard deviations {0} {1}\n".format(top, bottom))

    print('Calculating controls completed')
