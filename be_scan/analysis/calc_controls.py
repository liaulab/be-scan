"""
Author: Calvin XiaoYang Hu
Date: 231128

{Description: Calculate the control stats from the generated dataset}
"""

from pathlib import Path
import pandas as pd
from be_scan.plot._annotating_ import *

def calc_controls(annotated_lib, comparisons, 
                  neg_ctrl_col, neg_ctrl_conditions,

    out_dir='', out_file='stats.txt', 
    ):

    """[Summary]
    Calculate and output the normalization controls of the dataframe conditions. 

    Parameters
    ----------
    annotated_lib : str or path
        String or path to the csv file containing the values for comparison.
        The column headers must match the sample names in comparisons
    comparisons : list of str
        a list of columns of the annotated_lib for which to calculate negative controls
    neg_ctrl_col : str
        column of annotated_lib which correspond to normalization control
    neg_ctrl_conditions : list of str
        names of categories of neg_ctrl_col to normalize dataframe
    
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'conditions.csv'
        Name of output dataframe with guides and counts. 
    save : bool, default True
        Whether or not to save the resulting dataframe
    """

    path = Path.cwd()
    df_conds = pd.read_csv(annotated_lib)

    # normalize data to intergenic controls if neg_ctrl is provided
    # calculate negative control stats
    
    _, list_negctrlstats, _ = calc_neg_ctrls(df_conds, 
                                             comparisons, 
                                             neg_ctrl_col, 
                                             neg_ctrl_conditions)
    f = open("demofile2.txt", "a")
    f.write("Now the file has more content!")
    f.close()

    # write and export files
    outpath = path / out_dir
    Path.mkdir(outpath, exist_ok=True)

    f = open(outpath / out_file, "w")
    for comp, mean, stdev, top, bottom in list_negctrlstats: 
        f.write("For comparison {0}".format(comp))
        f.write("Mean is {0}".format(mean))
        f.write("Standard deviation is {0}".format(stdev))
        f.write("Mean +- 2 standard deviations {0} {1}".format(top, bottom))
    f.close()

    print('calc_controls outputed to', str(outpath / out_file))
    print('Calculating controls completed')
