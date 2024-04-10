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
import matplotlib as mpl

from be_scan.plot._annotating_ import *

def boxplot(df_filepath, comparisons, # each comparison is a plot
               plot_column, plot_conditions, # each plot condition is a box in a plot

    # filter out unwanted quantitative params
    filter_val=False, val_cols=['CtoT_muttypes'], val_min=0, 
    # filter out unwanted categorical params
    filter_params=False, 
    params_cols=['CtoT_muttype'], 
    params_conditions=[['Missense', 'Silent', 'Mixed', 'Nonsense']], 
    # normalization params
    neg_ctrl=False, neg_ctrl_col='CtoT_win_overlap', neg_ctrl_conditions=['Intron', 'Control'], 
    
    # boxplot labels
    xlab='', ylab='Log2 Fold Change', col_label='subavg', 
    # output params
    savefig=True, out_name='boxes', out_type='png', out_directory='', show=True, 

    subplots_kws={'figsize':(5,4)}, 
    boxplot_kws = {'saturation':1, 'fliersize':4, 'width':0.4, 
                   'flierprops':{'marker':'o', 'mec':'black', 'lw':1, 'alpha':0.8}},
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

    filter_val : bool, optional, defaults to False
        whether or not to exclude a subset of data from plotting by a minimum value
        default purpose is to filter out all intron and exon/intron guides
    val_cols : list of str, optional, 
        defaults to ['CtoT_muttypes']
        names of columns to filter dataframe for plotting
    val_min : int, optional, defaults to 0
        the minimum value by which to filter rows by val_cols

    filter_params : bool, optional, defaults to False
        whether or not to exclude a subset of data from plotting by categorical params
        default purpose is to filter to keep only Missense, Nonsense, Silent, Mixed guides
    params_cols : list of str, optional, 
        defaults to ['CtoT_muttype']
        names of column to filter dataframe for plotting
    params_conditions : list of lists of str, optional, 
        defaults to [['Missense', 'Silent', 'Mixed', 'Nonsense']]
        names of categories of filter_col to filter dataframe
        
    neg_ctrl : bool, optional, defaults to False
        whether or not to calulate negative control for normalization and line drawing
    neg_ctrl_col : str, optional, defaults to 'CtoT_win_overlap'
        column of .csv which correspond to normalization control
    neg_ctrl_conditions : list of str, optional, defaults to ['Intron', 'Control']
        name of categories of neg_ctrl_col to normalize dataframe

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
    show : bool, optional, defaults to True
        whether or not to show the plot

    subplots_kws: dict, optional, defaults to 
        {'figsize':(5,4)}
        input params for plt.subplots() 
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
    boxplot_kws: dict, optional, defaults to 
        {'saturation':1, 'fliersize':4, 'width':0.4, 
        'flierprops':{'marker':'o', 'mec':'black', 'lw':1, 'alpha':0.8}}
        input params for sns.boxplot() 
        https://seaborn.pydata.org/generated/seaborn.boxplot.html
    axhline_kws: dict, optional, defaults to 
        {'color':'k', 'ls':'--', 'lw':1}
        input params for plt.axhline() 
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axhline.html

    Returns
    ----------
    None
    """

    # style
    mpl.rcParams.update({'font.size': 10})

    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)

    # check conflicting params and output for user
    if filter_val: 
        assert isinstance(val_min, float), "check param: val_min"
        assert isinstance(val_cols, list) and len(val_cols) > 0, "check param: val_cols"
        for vc in val_cols: 
            assert vc in df_data.columns.tolist(), "check param: val_cols"
    if filter_params: 
        assert isinstance(params_cols, list), "check param: params_cols"
        assert isinstance(params_conditions, list), "check param: params_conditions"
        for pc in params_cols: 
            assert pc in df_data.columns.tolist(), "check param: val_cols"
    if neg_ctrl: 
        assert isinstance(neg_ctrl_col, str), "check param: params_cols"
        assert isinstance(neg_ctrl_conditions, list), "check param: params_conditions"
        assert neg_ctrl_col in df_data.columns.tolist(), "check param: val_cols"
    assert isinstance(subplots_kws, dict), "check param: subplots_kws"
    assert isinstance(boxplot_kws, dict), "check param: boxplot_kws"
    assert isinstance(axhline_kws, dict), "check param: axhline_kws"

    # normalize data to intergenic controls if neg_ctrl is provided
    if neg_ctrl: 
        # calculate negative control stats
        _, list_negctrlstats, avg_dict = calc_neg_ctrls(df_data, comparisons, 
                                                        neg_ctrl_col, neg_ctrl_conditions)
        # calculate normalized log_fc scores for each comp condition
        df_data = norm_to_intergenic_ctrls(df_data, comparisons, avg_dict, col_label)

    # filter out subset of data if filters are provided
    if filter_params: 
        for col, conds in zip(params_cols, params_conditions): 
            df_data = df_data.loc[df_data[col].isin(conds)]
    if filter_val: 
        for col in val_cols: 
            df_data = df_data[df_data[col] > val_min]
    df_data = df_data.loc[df_data[plot_column].isin(plot_conditions)].copy()
    
    for comp in comparisons:
        # make plot for every comparison
        _, ax = plt.subplots(**subplots_kws)
        y = comp+'_'+col_label if neg_ctrl else comp
        sns.boxplot(data=df_data, ax=ax, x=plot_column, y=y, 
                    **boxplot_kws)
        plt.setp(ax.artists, edgecolor='black')
        plt.setp(ax.lines, color='black')

        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        if neg_ctrl and list_negctrlstats != None:
            tup_comp_stdev = [tup for tup in list_negctrlstats if tup[0] == comp][0][2]
            plt.axhline(y=2*tup_comp_stdev, **axhline_kws) # top baseline
            plt.axhline(y=-2*tup_comp_stdev, **axhline_kws) # bottom baseline
        
        # Adjust x and y axis limits
        plt.ylim(np.floor(df_data[y].min()),
                 np.ceil(df_data[y].max()))
        # Set/adjust labels
        plt.title(comp)
        plt.ylabel(ylab) ; plt.xlabel(xlab)
        plt.xticks(rotation=45, horizontalalignment='right')
        # Adjust dimensions
        plt.tight_layout()

        # Save to pdf
        path = Path.cwd()
        if savefig: 
            outpath = path / out_directory
            out = out_name + comp + '.' + out_type
            plt.savefig(outpath / out, format=out_type)
        if show: 
            plt.show()
        plt.close()
