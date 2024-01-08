"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function performs normalization and then plots the data for each condition to reveal enriched guides}
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
from be_scan.plot._annotating_ import *

def plot_scatterplot(df_filepath, # dataframe
                     x_column, y_column, 
                     hue_column, # df params
                     comparisons, 
                     neg_ctrl_col, neg_ctrl_category,
                     xmin=None, xmax=None,
                     xlab='Amino Acid Position', ylab='sgRNA Score', # scatterplot labels
                     out_name='scatterplot', out_type='pdf', out_directory='', # output params
                     savefig=True,
                     scatterplot_kws={'alpha':0.8, 'linewidth':1.0, 
                                         'edgecolor':'black', 's':25},
                     subplots_kws={'figsize':(8,4)},
                     axhline_kws={'color':'k', 'ls':'--', 'lw':1},
                     autoannotate=False,
                     autoannotate_label=None, autoannotate_top=None, autoannotate_cutoff=None,
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
    y_column : str, required
        column of .csv, typically the normalized log_fc score
    hue_column : str, required
        column of .csv which correspond to coloring of points
    comparisons : list of str, required
        list of comparisons that correspond to columns of data
    neg_ctrl_col : str, required
        column of .csv which correspond to normalization control
    neg_ctrl_category : str, required
        categorical variable of neg_ctrl_col

    xlab : str, optional, defaults to 'Amino Acid Position'
        name of x-axis label
    ylab : str, optional, defaults to 'sgRNA Score'
        name of y-axis label
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
        ath to output directory
    
    autoannotate : bool, False
        Whether or not to autoannotate points
    autoannotate_label : string, None
        The column of the label in the dataframe
    autoannotate_top : int, None
        The top n scoring points will be labeled
    autoannotate_cutoff : float, None
        The absolute value cutoff for which points will be labeled

    Returns
    ----------
    None
    """

    in_dataframe = pd.read_csv(df_filepath)
    # calculate negative control stats
    _, list_negctrlstats, avg_dict = calc_negative_controls(in_dataframe, comparisons, neg_ctrl_col, neg_ctrl_category)
    # calculate normalized log_fc scores for each comp condition
    df_logfc = norm_to_intergenic_ctrls(in_dataframe, comparisons, avg_dict, y_column)

    hue_order_s = list_muttypes[:3]
    palette_s = color_list[:3]

    df_filtered = df_logfc.loc[df_logfc[hue_column].isin(list_muttypes[:3])]
    
    for comp in comparisons:
        # Make plots
        fig, ax = plt.subplots(**subplots_kws)
        # name of y values is comp + y_column
        y = comp+'_'+y_column 
        # scatterplot
        sns.scatterplot(ax=ax,
                        data=df_filtered, # dataframe of edits in window, normalized
                        x=x_column, y=y, # x and y columns plotted against each other
                        hue=hue_column, # color scatterplot according to column
                        hue_order=hue_order_s, palette=palette_s, 
                        **scatterplot_kws
                        )
        
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        tup_comp_stdev = [tup for tup in list_negctrlstats if tup[0] == comp][0][2]
        ax.axhline(y=2*tup_comp_stdev, **axhline_kws) # top baseline
        ax.axhline(y=-2*tup_comp_stdev, **axhline_kws) # bottom baseline

        # autoannotate the scoring hits
        if autoannotate: 
            if autoannotate_cutoff: 
                for _, row in df_filtered[df_filtered[y].abs() > autoannotate_cutoff].iterrows():
                    if row[autoannotate_label] != "nan": 
                        plt.text(row[x_column], row[y], row[autoannotate_label], ha='right', va='bottom', fontsize=6)
            elif autoannotate_top: 
                print('Second')
                for _, row in df_filtered.assign(yabs=df_filtered[y].abs()).sort_values(by='yabs', ascending=False).head(autoannotate_top).iterrows():
                    plt.text(row[x_column], row[y], row[autoannotate_label], ha='right', va='bottom', fontsize=6)
            else: 
                print("Please include a cutoff or top value")
        
        # Adjust x and y axis limits
        ax.set_xlim(df_filtered[x_column].min()-10, df_filtered[x_column].max()+10)
        bound = max(abs(np.floor(df_filtered[y].min())), 
                    abs(np.ceil(df_filtered[y].max())))
        ax.set_ylim(-1*bound, bound)
        # Set title and axes labels
        ax.set_title(comp)
        ax.set_xlim(xmin, xmax)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)

        # Adjust figsize
        plt.tight_layout()
        # Save to pdf
        if savefig:
            output_path = out_directory + out_name + comp + '.' + out_type
            plt.savefig(output_path, format=out_type)
        plt.show()
        plt.close()
