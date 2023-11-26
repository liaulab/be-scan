"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from be_scan.plot._annotating_ import list_muttypes, color_list

def plot_corr_scatterplot(df_filepath, 
                          condition1, condition2, 
                          hue_column, 
                          hue_order=list_muttypes, palette=color_list, 
                          xmin=None, xmax=None, ymin=None, ymax=None, 
                          xlab='cond1 score', ylab='cond2 score', 
                          out_directory='', out_name='correlation_scatterplot', out_type='pdf', 
                          alpha=0.8, linewidth=1, edgecolor='black', s=25,
                          figsize=(4.5, 4), 
                          savefig=True,
                          ):
    
    """[Summary]
    This function takes in a dataframe from count_reads, and plots
    a heatmap showing correlation between all given comparison conditions
    ...
    
    :param df_filepath: filepath to .csv data generated from count_reads
    :type df_filepath: str, required
    :param condition1: comparison condition 1, name of a column in .csv data
    :type condition1: str, required
    :param condition2: comparison condition 2, name of a column in .csv data
    :type condition2: str, required
    :param hue_column: the categorial data for each of the points, name of a column in .csv data
    :type hue_column: str, required

    :param hue_order: a list of categorial variables in hue_column
    :type hue_order: list of str, optional, defaults to list_muttypes in _annotating_.py a preset list of column names
    :param palette: a list of colors which correspond to hue_order
    :type palette: list of str, optional, defaults to color_list in _annotating_.py a preset list of colors from ColorBrewer2

    :param xmin: x-axis left bound
    :type xmin: float, optional, defaults to None
    :param xmax: x-axis right bound
    :type xmax: float, optional, defaults to None
    :param ymin: y-axis lower bound
    :type ymin: float, optional, defaults to None
    :param ymax: y-axis upper bound
    :type ymax: float, optional, defaults to None
    :param xlab: x-axis label
    :type xlab: str, optional, defaults to 'cond1 score'
    :param ylab: y-axis label
    :type ylab: str, optional, defaults to 'cond2 score'

    :param out_name: name of the output plot
    :type out_name: str, optional, defaults to 'scatterplot'
    :param out_type: type of the output plot
    :type out_type: str, optional, defaults to 'pdf'
    :param out_directory: directory path of the output plot
    :type out_directory: str, optional, defaults to ''

    :param alpha: transparency of scatterplot points
    :type alpha: float, optional, defaults to 0.8
    :param linewidth: linewidth of plot
    :type linewidth: float, optional, defaults to 1.0
    :param edgecolor: color of scatterplot edge lines
    :type edgecolor: str, optional, defaults to 'black'
    :param s: size of scatterplot points
    :type s: int, optional, defaults to 25
    :param dimensions: the figsize (length, width)
    :type dimensions: tuple of ints, optional, defaults to (8,4)
    :param savefig: option of saving figure to output or not
    :type figsize: boolean, optional, defaults to True
    ...
    
    :return: None
    :rtype: NoneType
    """

    df_data = pd.read_csv(df_filepath)
    df_filtered = df_data.loc[df_data[hue_column].isin(hue_order)]
    
    # Make plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(data=df_filtered, 
                    ax=ax, 
                    x=condition1, y=condition2, 
                    hue=hue_column, hue_order=hue_order, palette=palette, 
                    alpha=alpha, linewidth=linewidth, edgecolor=edgecolor, s=s
                    )
    
    # Adjust x and y axis limits
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    # Set labels
    plt.xlabel(xlab) # set x-axis label
    plt.ylabel(ylab) # set y-axis label
    # Adjust dimensions
    plt.tight_layout()
    plt.show()

    # Save to pdf
    out = out_directory + condition1 + condition2 + '_' + out_name + '.' + out_type
    if savefig: 
        plt.savefig(out, format='pdf')
    plt.close()

# python3 -m be_scan plot_corr_scatterplot -df '../../../Downloads/NZL10196_v9_comparisons.csv' -c1 'd3-neg' -c2 'd9-pos' -hue 'Mut_type'
