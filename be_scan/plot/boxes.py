"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function plots chosen guides by plot_column categories to show the distribution of guides}
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from be_scan.plot._annotating_ import norm_to_intergenic_ctrls, calc_negative_controls

def plot_boxes(df_filepath, 
               comparisons, # each comparison is a plot
               plot_column, plot_conditions, # each plot condition is a box in a plot

               neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # normalization params
               filter=False, filter_col='', filter_conditions=[], # filter out unwanted data params
               xlab='', ylab='Log2 Fold Change', col_label='subavg', # labels
               savefig=True, out_name='boxes', out_type='pdf', out_directory='', # output params

               subplots_kws={'figsize':(5,4)}, 
               boxplot_kws = {'saturation':1, 'fliersize':4, 'width':0.4, 
                              'flierprops':{'marker':'o', 'mec':'black', 'lw':1, 'alpha':0.8}
                              },
               axhline_kws = {'color':'k', 'ls':'--', 'lw':1},
               ):
    
    """[Summary]
    This function takes in a dataframe from count_reads, performs normalization, 
    and then plots chosen (ex control) guides by plot_column categories, 
    to show the distribution of categories of guides
    
    Parameters
    ----------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    comparisons : list of str, required
        list of comparisons that correspond to columns of data
    plot_column : str, required
        column of .csv, typically domain or mutation type
    plot_conditions : list of str, required
        category names of plot_column

    neg_ctrl : bool, optional, defaults to False
        whether or not to calulate negative control for normalization and line drawing
    neg_ctrl_col : str, optional, defaults to ''
        column of .csv which correspond to normalization control
    neg_ctrl_conditions : list of str, optional, defaults to []
        name of categories of neg_ctrl_col to normalize dataframe
    filter : bool, optional, defaults to False
        whether or not to exclude a subset of data from plotting
    filter_col : list of str, optional, defaults to []
        names of column to filter dataframe for plotting
    filter_conditions : list of lists of str, optional, defaults to []
        names of categories of filter_col to filter dataframe

    xlab : str, optional, defaults to ''
        name of x-axis label
    ylab : str, optional, defaults to 'Log2 Fold Change'
        name of y-axis label
    col_label : str, optional, defaults to 'subavg'
        a suffix label for point values adjusted by normalization

    savefig : boolean, optional, defaults to True
        option of saving figure to output or not
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
        path to output directory

    axhline_kws: dict, optional, defaults to {'color':'k', 'ls':'--', 'lw':1}
        input params for plt.axhline()
    boxplot_kws: dict, optional, defaults to {'saturation':1, 'fliersize':4, 'width':0.4, 
                                              'flierprops':{'marker':'o', 'mec':'black', 
                                              'lw':1, 'alpha':0.8}
        input params for sns.boxplot()
    subplots_kws: dict, optional, defaults to {'figsize':(4.5, 4)}
        input params for plt.subplots()

    Returns
    ----------
    None
    """
    
    df_filepath = Path(df_filepath)
    df_input = pd.read_csv(df_filepath)
    # normalize data to intergenic controls if neg_ctrl is provided
    if neg_ctrl: 
        # calculate negative control stats
        _, list_negctrlstats, avg_dict = calc_negative_controls(df_input, comparisons, neg_ctrl_col, neg_ctrl_conditions)
        # calculate normalized log_fc scores for each comp condition
        df_input = norm_to_intergenic_ctrls(df_input, comparisons, avg_dict, col_label)
    # filter out subset of data if filter is provided
    if filter: 
        for col, conds in zip(filter_col, filter_conditions): 
            df_input = df_input.loc[df_input[col].isin(conds)]
    df_input = df_input.loc[df_input[plot_column].isin(plot_conditions)].copy()
    
    for comp in comparisons:
        # make plot for every comparison
        _, ax = plt.subplots(**subplots_kws)
        y = comp+'_'+col_label if neg_ctrl else comp
        sns.boxplot(data=df_input, ax=ax, x=plot_column, y=y, 
                    **boxplot_kws
                    )
        plt.setp(ax.artists, edgecolor='black')
        plt.setp(ax.lines, color='black')

        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        if neg_ctrl and list_negctrlstats != None:
            tup_comp_stdev = [tup for tup in list_negctrlstats if tup[0] == comp][0][2]
            plt.axhline(y=2*tup_comp_stdev, **axhline_kws) # top baseline
            plt.axhline(y=-2*tup_comp_stdev, **axhline_kws) # bottom baseline
        
        # Adjust x and y axis limits
        plt.ylim(np.floor(df_input[y].min()),
                 np.ceil(df_input[y].max()))
        # Set/adjust labels
        plt.title(comp) # Set plot title
        plt.ylabel(ylab) # Set y-axis label
        plt.xlabel(xlab) # Remove x-axis label
        plt.xticks(rotation=45, horizontalalignment='right')
        
        # Adjust dimensions
        plt.tight_layout()
        # Save to pdf
        path = Path.cwd()
        if savefig: 
            outpath = path / out_directory
            out = out_name + comp + '.' + out_type
            plt.savefig(outpath / out, format=out_type)
        plt.show()
        plt.close()
