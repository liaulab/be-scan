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

def plot_scatterplot(df_filepath, # dataframe
                     comparisons, # each comparison is a plot, and also the y axis
                     x_column, # the x axis values
                     
                     hue=False, hue_column='', palette=color_list, # color params
                     neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # normalization params
                     filter=False, filter_col=[], filter_conditions=[], # filter out unwanted data params
                     autoannot=False, autoannot_label=None, autoannot_top=None, autoannot_cutoff=None, # autoannotate outliers
                     xlab='Amino Acid Position', ylab='sgRNA Score', col_label='subavg', # scatterplot labels
                     savefig=True, out_name='scatterplot', out_type='pdf', out_directory='', # output params

                     xlim_kws={'xmin':None, 'xmax':None}, ylim_kws={'ymin':None, 'ymax':None},
                     scatterplot_kws={'alpha':0.8, 'linewidth':1.0, 
                                         'edgecolor':'black', 's':25},
                     subplots_kws={'figsize':(8,4)},
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
        
    hue : bool, optional, defaults to False
        whether or not to color points by a column value
    hue_column : str, optional, defaults to ''
        column of .csv which correspond to coloring of points
    palette : list, optional, defaults to color_list in be_scan.plot._annotating_
        list of color codes

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
        ath to output directory

    xlim_kws: dict, optional, defaults to {'xmin':None, 'xmax':None}
        input params for ax.set_xlim()
    ylim_kws: dict, optional, defaults to {'ymin':None, 'ymax':None}
        input params for ax.set_ylim()
    scatterplot_kws: dict, optional, defaults to {'alpha':0.8, 'linewidth':1, 'edgecolor':'black', 's':25}
        input params for sns.scatterplot()
    subplots_kws: dict, optional, defaults to {'figsize':(4.5, 4)}
        input params for plt.subplots()

    Returns
    ----------
    None
    """

    df_input = pd.read_csv(df_filepath)
    # normalize data to intergenic controls if neg_ctrl is provided
    if neg_ctrl: 
        # calculate negative control stats
        _, list_negctrlstats, avg_dict = calc_negative_controls(df_input, comparisons, neg_ctrl_col, neg_ctrl_conditions)
       # calculate normalized log_fc scores for each comp condition
        df_input = norm_to_intergenic_ctrls(df_input, comparisons, avg_dict, col_label)
    if filter: 
        for col, conds in zip(filter_col, filter_conditions): 
            df_input = df_input.loc[df_input[col].isin(conds)]
    
    for comp in comparisons:
        # Make plots
        _, ax = plt.subplots(**subplots_kws)
        y = comp+'_'+col_label if neg_ctrl else comp
        baseline_params = {'data':df_input, 'ax':ax, 'x':x_column, 'y':y}
        if hue: 
            sns.scatterplot(**baseline_params, 
                            hue=hue_column, palette=palette, 
                            **scatterplot_kws
                            )
        else: 
            sns.scatterplot(**baseline_params, 
                            **scatterplot_kws
                            )
        
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        if neg_ctrl and list_negctrlstats != None:
            tup_comp_stdev = [tup for tup in list_negctrlstats if tup[0] == comp][0][2]
            plt.axhline(y=2*tup_comp_stdev, **axhline_kws) # top baseline
            plt.axhline(y=-2*tup_comp_stdev, **axhline_kws) # bottom baseline

        # autoannot the scoring hits
        if autoannot: 
            if autoannot_cutoff: # annotate according to a cutoff where all values above are annotated
                for _, row in df_input[df_input[y].abs() > autoannot_cutoff].iterrows():
                    if row[autoannot_label] != "nan": 
                        plt.text(row[x_column], row[y], row[autoannot_label], ha='right', va='bottom', fontsize=6)
            elif autoannot_top: # annotate according to the top n scoring
                for _, row in df_input.assign(yabs=df_input[y].abs()).sort_values(by='yabs', ascending=False).head(autoannot_top).iterrows():
                    plt.text(row[x_column], row[y], row[autoannot_label], ha='right', va='bottom', fontsize=6)
            else: 
                print("Please include a cutoff or top value")
        
        # Adjust x and y axis limits
        ax.set_xlim(df_input[x_column].min()-10, df_input[x_column].max()+10)
        bound = max(abs(np.floor(df_input[y].min())), 
                    abs(np.ceil(df_input[y].max())))
        ax.set_ylim(-1*bound, bound)
        # Set title and axes labels
        ax.set_title(comp)
        ax.set_xlim(**xlim_kws)
        ax.set_ylim(**ylim_kws)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)

        # Adjust figsize
        plt.tight_layout()
        # Save to pdf
        if savefig:
            out = out_name + comp + '.' + out_type
            output_path = Path(out_directory)
            plt.savefig(output_path / out, format=out_type)
        plt.show()
        plt.close()
