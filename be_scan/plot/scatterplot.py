"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from ._annotating_ import color_list, list_muttypes
from ._annotating_ import norm_to_intergenic_ctrls, calc_negative_controls

def plot_scatterplot(df_filepath, # dataframe
                     x_column, y_column, 
                     hue_column, # df params
                     comparisons, 
                     neg_ctrl_col, neg_ctrl_category,
                     window=None,
                     xlab='Amino Acid Position', ylab='sgRNA Score', # scatterplot labels
                     out_name='scatterplot', out_type='pdf', out_directory='', # output params
                     alpha=0.8, linewidth=1.0, edgecolor='black', s=25, # scatterplot params
                     figsize=(8,4)
                     ):
    
    """[Summary]
    This function takes in a dataframe from count_reads, performs normalization, 
    and then plots the data for each condition to reveal which guides are enriched
    ...

    :param df_filepath: filepath to .csv data generated from count_reads
    :type df_filepath: str, required
    :param x_column: column of .csv for x axis, typically amino acid position
    :type x_column: str, required
    :param y_column: column of .csv for y axis, typically the normalized log_fc change score
    :type y_column: str, required
    :param hue_column: column of .csv which correspond to coloring of scatterplot points
    :type hue_column: str, required
    :param comparisons: list of comparisons that correspond to columns of .csv data
    :type comparisons: list of str, required
    :param neg_ctrl_col: column of .csv which correspond to normalization control
    :type neg_ctrl_col: str, required
    :param neg_ctrl_category: categorical variable of neg_ctrl_col .csv which correspond to normalization control
    :type neg_ctrl_category: str, required

    :param window: inclusive window of which amino acid positions are shown in the plot
    :type window: tuple of ints, optional, defaults to None
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

    df_filtered = df_logfc.loc[(df_logfc[x_column]>=window[0]) & (df_logfc[x_column]<=window[1])]
    df_filtered = df_filtered.loc[df_filtered[hue_column].isin(list_muttypes[:3])]

    # output pdf information
    output_path = out_directory + out_name + '.' + out_type
    figpdf = PdfPages(output_path)
    
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
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)

        # Adjust figsize
        plt.tight_layout()
        # Save to pdf
        plt.savefig(figpdf, format='pdf')
        plt.close()

    figpdf.close()

# python3 -m be_scan plot_scatterplot -df '../../../Downloads/NZL10196_v9_comparisons.csv' 
#         -x 'Edit_site_3A1' -y 'log2_fc' -c 'd3-pos' -hue 'Mut_type' -pt 'comparison' 
#         -ncol 'Gene' -ncat 'NON-GENE' -win 224 912 -c 'd3-pos' 'd3-neg' 'd6-pos' 'd6-neg' 'd9-pos' 'd9-neg'
