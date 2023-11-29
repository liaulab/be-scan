"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function plots a heatmap showing correlation between all conditions}
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def plot_corr_heatmap(df_filepath, 
                      comparisons, 
                      corr_type='spearman', 
                      xlab='', ylab='', title='Spearman Correlation Heatmap', # output related params
                      out_directory='', out_name='correlation_heatmap', out_type='pdf', 
                      center=0, linewidth=0.5, line_pos=[], cmap='coolwarm', # plot related params
                      savefig=True,
                     ): 
    
    """[Summary]
    This function takes in a dataframe from count_reads, and plots
    a scatterplot showing correlation between two given conditions
    ...
    
    :param df_filepath: filepath to .csv data generated from count_reads
    :type df_filepath: str, required
    :param comparisons: list of comparisons that correspond to columns of .csv data
    :type comparisons: list of str, required

    :param corr_type: type of correlation calculation
    :type corr_type: str, optional, defaults to 'spearman'
    :param xlab: name of the x-axis label
    :type xlab: str, optional, defaults to ''
    :param ylab: name of the y-axis label
    :type ylab: str, optional, defaults to ''
    :param title: name of title label
    :type title: str, optional, defaults to 'Spearman Correlation Heatmap'

    :param out_directory: path to output directory
    :type out_directory: str, optional, defaults to ''
    :param out_name: name of figure output
    :type out_name: str, optional, defaults to 'scatterplot'
    :param out_type: file type of figure output
    :type out_type: str, optional, defaults to 'pdf'

    :param center: value for centering the colormap
    :type center: float, optional, defaults to 0
    :param linewidth: linewidth of lines dividing cells
    :type linewidth: float, optional, defaults to 0.5
    :param line_pos: positions to put larger lines for dividing the cells
    :type line_pos: list of int, optional, defaults to []
    :param cmap: matplotlib colormap
    :type cmap: str, optional, defaults to 'coolwarm'
    :param savefig: option of saving figure to output or not
    :type figsize: boolean, optional, defaults to True
    ...
    
    :return: None
    :rtype: NoneType
    """

    # Plotting parameters and variables
    sns.set_style('ticks')

    # Compute correlation matrix
    df = pd.read_csv(df_filepath)
    df_comp = df[comparisons].copy()
    df_corr = df_comp.corr(method=corr_type)

    # Set up the matplotlib figure
    _, ax = plt.subplots(figsize=(4,4))
    ax = sns.heatmap(df_corr, 
                     cmap=cmap, center=center,
                     square=True, linewidths=linewidth, cbar_kws={"shrink": 0.5}
                     )
    
    # show frame
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    
    # rotate axis labels
    if line_pos: 
        for pos in line_pos:
            plt.hlines(pos, *ax.get_xlim())
            plt.vlines(pos, *ax.get_ylim())
    # adjustments and labels
    plt.title(title)
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.yticks(rotation=0, horizontalalignment='right')
    plt.ylabel(xlab)
    plt.xlabel(ylab)
    plt.tight_layout()

    # save pdf and close everything
    if savefig: 
        output_path = out_directory + out_name + '.' + out_type
        plt.savefig(output_path, format=out_type)
    plt.show()
    plt.close()
