"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function performs normalization and then plots the data for each condition to reveal enriched guides}
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

import pandas as pd
from be_scan.plot._annotating_ import *

mpl.rcParams.update({'font.size': 10})
cm = 1/2.54
fig,ax = plt.subplots()
fig.set_size_inches(9*cm,7*cm)

def scatterplot(df_filepath, # dataframe
                comparisons, # each comparison is a plot, and also the y axis
                x_column, # the x axis values
                     
    filter_val=False, val_cols=[], val_min=None, # filter out unwanted quantitative params
    filter_params=False, params_cols=[], params_conditions=[], # filter out unwanted categorical params
    include_hue=False, hue_col='Mut_type', hue_order=list_muttypes, palette=color_list, # color params
    neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # normalization params

    autoannot=False, autoannot_label=None, autoannot_top=None, autoannot_cutoff=None, # autoannotate outliers
    xlab='Amino Acid Position', ylab='sgRNA Score', col_label='subavg', # scatterplot labels
    savefig=True, out_name='scatterplot', out_type='png', out_directory='', show=True, # output params

    xlim_kws={'xmin':None, 'xmax':None}, ylim_kws={'ymin':None, 'ymax':None},
    scatterplot_kws={'alpha':0.8, 'linewidth':1.0, 
                        'edgecolor':'black', 's':25},
    subplots_kws={'figsize':(10,4)},
    axhline_kws={'color':'k', 'ls':'--', 'lw':1},
    ):
    
    """[Summary]
    This function takes in a dataframe from count_reads, performs normalization, 
    and then plots the data for each condition to reveal which guides are enriched

    Parameters
    ----------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    x_column : str, required
        column of .csv, typically amino acid position
    comparisons : list of str, required
        list of comparisons that correspond to columns of data

    filter_val : bool, optional, defaults to False
        whether or not to exclude a subset of data from plotting by a minimum value
    val_cols : list of str, optional, defaults to []
        names of columns to filter dataframe for plotting
    val_min : list of str, optional, defaults to None
        the minimum value by which to filter rows by val_cols
    filter_params : bool, optional, defaults to False
        whether or not to exclude a subset of data from plotting by categorical params
    params_cols : list of str, optional, defaults to []
        names of column to filter dataframe for plotting
    params_conditions : list of lists of str, optional, defaults to []
        names of categories of filter_col to filter dataframe
    include_hue: bool, optional, default to False
        whether or not to color points by a variable, 
        will also restrict points plotted to only the hue_order values listed
    hue_col: str, optional, defaults to 'Mut_type'
        the categorial dimension of the data, name of .csv data column
    hue_order: list of str, optional, defaults to a preset list of column names
        a list of categorial variables in hue_col
    palette: list of str, optional, defaults to a preset list of colors from ColorBrewer2
        a list of colors which correspond to hue_order
    neg_ctrl : bool, optional, defaults to False
        whether or not to calulate negative control for normalization and line drawing
    neg_ctrl_col : str, optional, defaults to ''
        column of .csv which correspond to normalization control
    neg_ctrl_conditions : list of str, optional, defaults to []
        name of categories of neg_ctrl_col to normalize dataframe

    autoannot : bool, False
        Whether or not to autoannot points
    autoannot_label : string, None
        The column of the label in the dataframe
    autoannot_top : int, None
        The top n scoring points will be labeled
    autoannot_cutoff : float, None
        The absolute value cutoff for which points will be labeled
    
    xlab : str, optional, defaults to 'Amino Acid Position'
        name of x-axis label
    ylab : str, optional, defaults to 'sgRNA Score'
        name of y-axis label
    col_label : str, optional, defaults to 'subavg'
        a suffix label for point values adjusted by normalization
    savefig: bool, optional, defaults to True
        whether or not to save the figure
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
        path to output directory
    show : bool, optional, defaults to True
        whether or not to show the plot

    xlim_kws: dict, optional, defaults to 
        {'xmin':None, 'xmax':None}
        input params for ax.set_xlim() 
        https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xlim.html
    ylim_kws: dict, optional, defaults to 
        {'ymin':None, 'ymax':None}
        input params for ax.set_ylim() 
        https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_ylim.html
    scatterplot_kws: dict, optional, defaults to 
        {'alpha':0.8, 'linewidth':1.0, 'edgecolor':'black', 's':25}
        input params for sns.scatterplot() 
        https://seaborn.pydata.org/generated/seaborn.scatterplot.html
    subplots_kws: dict, optional, defaults to 
        {'figsize':(4,4)}
        input params for plt.subplots() 
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
    axhline_kws: dict, optional, defaults to 
        {'color':'k', 'ls':'--', 'lw':1}
        input params for plt.axhline() 
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axhline.html

    Returns
    ----------
    None
    """

    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)

    # check conflicting params and output for user
    if filter_val: 
        assert isinstance(val_min, float), "check param: val_min"
        assert isinstance(val_cols, list) and len(val_cols) > 0, "check param: val_cols"
    if filter_params: 
        assert isinstance(params_cols, list), "check param: params_cols"
        assert isinstance(params_conditions, list), "check param: params_conditions"
    if include_hue: 
        assert hue_col in df_data.columns.tolist(), "check param: hue_col"
        assert isinstance(hue_order, list) and len(hue_order) > 0, "check param: hue_order"
        assert isinstance(palette, list) and len(palette) > 0, "check param: palette"
    if neg_ctrl: 
        assert isinstance(neg_ctrl_col, str), "check param: params_cols"
        assert isinstance(neg_ctrl_conditions, list), "check param: params_conditions"
        assert neg_ctrl_col in df_data.columns.tolist(), "check param: val_cols"
    if autoannot: 
        assert isinstance(autoannot_label, str), "check param: autoannot_label"
        assert autoannot_label in df_data.columns.tolist(), "check param: autoannot_label"
        autoannot_bool = isinstance(autoannot_top, int) or isinstance(autoannot_cutoff, float)
        assert autoannot_bool, "check param: autoannot_top autoannot_cutoff"
    assert isinstance(xlim_kws, dict), "check param: xlim_kws"
    assert isinstance(ylim_kws, dict), "check param: ylim_kws"
    assert isinstance(scatterplot_kws, dict), "check param: scatterplot_kws"
    assert isinstance(subplots_kws, dict), "check param: subplots_kws"
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

    for comp in comparisons:
        # make plots
        _, ax = plt.subplots(**subplots_kws)
        y = comp+'_'+col_label if neg_ctrl else comp
        baseline_params = {'data':df_data, 'ax':ax, 'x':x_column, 'y':y}
        if include_hue: 
            baseline_params = {**baseline_params, 'hue':hue_col, 'palette':palette}
        sns.scatterplot(**baseline_params, **scatterplot_kws)
        
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        if neg_ctrl and list_negctrlstats != None:
            tup_comp_stdev = [tup for tup in list_negctrlstats if tup[0] == comp][0][2]
            plt.axhline(y=2*tup_comp_stdev, **axhline_kws) # top baseline
            plt.axhline(y=-2*tup_comp_stdev, **axhline_kws) # bottom baseline

        # autoannot the screen hits
        if autoannot: 
            if autoannot_cutoff: # annotate according to a cutoff where all values above are annotated
                for _, row in df_data[df_data[y].abs() > autoannot_cutoff].iterrows():
                    if row[autoannot_label] != "nan": 
                        plt.text(row[x_column], row[y], row[autoannot_label], 
                                 ha='right', va='bottom', fontsize=6)
            elif autoannot_top: # annotate according to the top n scoring
                sorted_df = df_data.assign(yabs=df_data[y].abs()).sort_values(by='yabs', ascending=False)
                toprows_df = sorted_df.head(autoannot_top)
                for _, row in toprows_df.iterrows(): ###
                    plt.text(row[x_column], row[y], row[autoannot_label], 
                             ha='right', va='bottom', fontsize=6)
            else: 
                print("Please include a cutoff or top value")
        
        # Adjust x and y axis limits
        ax.set_xlim(df_data[x_column].min()-10, df_data[x_column].max()+10)
        bound = max(abs(np.floor(df_data[y].min())), 
                    abs(np.ceil(df_data[y].max())))
        ax.set_ylim(-1*bound, bound)
        # Set title and axes labels
        ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        ax.set_title(comp)
        ax.set_xlim(**xlim_kws) ; ax.set_ylim(**ylim_kws) ### set_xlim repeated
        ax.set_xlabel(xlab) ; ax.set_ylabel(ylab)
        # Adjust figsize
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
