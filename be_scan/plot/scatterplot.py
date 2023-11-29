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
from be_scan.plot._annotating_ import color_list, list_muttypes
from be_scan.plot._annotating_ import norm_to_intergenic_ctrls, calc_negative_controls

def plot_scatterplot(df_filepath, # dataframe
                     x_column, y_column, 
                     hue_column, # df params
                     comparisons, 
                     neg_ctrl_col, neg_ctrl_category,
                     xmin=None, xmax=None,
                     xlab='Amino Acid Position', ylab='sgRNA Score', # scatterplot labels
                     out_name='scatterplot', out_type='pdf', out_directory='', # output params
                     savefig=True,
                     alpha=0.8, linewidth=1.0, edgecolor='black', s=25, # scatterplot params
                     figsize=(8,4), 
                     ):
    
    """[Summary]
    This function takes in a dataframe from count_reads, performs normalization, 
    and then plots the data for each condition to reveal which guides are enriched
    ...

    :param df_filepath: filepath to .csv data generated from count_reads
    :type df_filepath: str, required
    :param x_column: column of .csv, typically amino acid position
    :type x_column: str, required
    :param y_column: column of .csv, typically the normalized log_fc score
    :type y_column: str, required
    :param hue_column: column of .csv which correspond to coloring of points
    :type hue_column: str, required
    :param comparisons: list of comparisons that correspond to columns of data
    :type comparisons: list of str, required
    :param neg_ctrl_col: column of .csv which correspond to normalization control
    :type neg_ctrl_col: str, required
    :param neg_ctrl_category: categorical variable of neg_ctrl_col
    :type neg_ctrl_category: str, required

    :param xlab: name of x-axis label
    :type xlab: str, optional, defaults to 'Amino Acid Position'
    :param ylab: name of y-axis label
    :type ylab: str, optional, defaults to 'sgRNA Score'
    :param out_name: name of figure output
    :type out_name: str, optional, defaults to 'scatterplot'
    :param out_type: file type of figure output
    :type out_type: str, optional, defaults to 'pdf'
    :param out_directory: path to output directory
    :type out_directory: str, optional, defaults to ''
    
    :param alpha: transparency of scatterplot points
    :type alpha: float, optional, defaults to 0.8
    :param linewidth: linewidth of plot
    :type linewidth: float, optional, defaults to 1.0
    :param edgecolor: color of scatterplot edge lines
    :type edgecolor: str, optional, defaults to 'black'
    :param s: size of scatterplot points
    :type s: int, optional, defaults to 25
    :param figsize: the figsize (length, width)
    :type figsize: tuple of ints, optional, defaults to (8,4)
    :param savefig: option of saving figure to output or not
    :type figsize: boolean, optional, defaults to True
    ...

    :return: None
    :rtype: NoneType
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
        fig, ax = plt.subplots(figsize=figsize)

        # name of y values is comp + y_column
        y = comp+'_'+y_column 
        # scatterplot
        sns.scatterplot(ax=ax,
                        data=df_filtered, # dataframe of edits in window, normalized
                        x=x_column, y=y, # x and y columns plotted against each other
                        hue=hue_column, # color scatterplot according to column
                        hue_order=hue_order_s, palette=palette_s, 
                        alpha=alpha, linewidth=linewidth, edgecolor=edgecolor, s=s)
        
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        tup_comp_stdev = [tup for tup in list_negctrlstats if tup[0] == comp][0][2]
        ax.axhline(y=2*tup_comp_stdev, color=edgecolor, ls='--', lw=linewidth) # top baseline
        ax.axhline(y=-2*tup_comp_stdev, color=edgecolor, ls='--', lw=linewidth) # bottom baseline
        
        # Adjust x and y axis limits
        ax.set_xlim(df_filtered[x_column].min()-10, df_filtered[x_column].max()+10)
        bound = max(abs(np.floor(df_filtered[y].min())), 
                    abs(np.ceil(df_filtered[y].max())))
        ax.set_ylim(-1*bound, bound)
        # Set title and axes labels
        ax.set_title(comp)
        ax.set_xlim(xmin,xmax)
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
