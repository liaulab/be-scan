"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function plots a scatterplot showing correlation between two given conditions}
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
                          scatterplot_kws={'alpha':0.8, 'linewidth':1, 
                                           'edgecolor':'black', 's':25}, 
                          subplots_kws = {'figsize':(4.5, 4)},
                          savefig=True,
                          ):
    
    """[Summary]
    This function takes in a dataframe from count_reads, and plots
    a heatmap showing correlation between all given comparison conditions
    
    Parameters
    ----------
    df_filepath: str, required
        filepath to .csv data generated from count_reads
    condition1: str, required
        comparison condition 1, name of .csv data column
    condition2: str, required
        comparison condition 2, name of .csv data column
    hue_column: str, required
        the categorial dimension of the data, name of .csv data column

    hue_order: list of str, optional, defaults to list_muttypes in _annotating_.py a preset list of column names
        a list of categorial variables in hue_column
    palette: list of str, optional, defaults to color_list in _annotating_.py a preset list of colors from ColorBrewer2
        a list of colors which correspond to hue_order

    xmin: float, optional, defaults to None
        x-axis left bound
    xmax: float, optional, defaults to None
        x-axis right bound
    ymin: float, optional, defaults to None
        y-axis lower bound
    ymax: float, optional, defaults to None
        y-axis upper bound
    xlab: str, optional, defaults to 'cond1 score'
        x-axis label
    ylab: str, optional, defaults to 'cond2 score'
        y-axis label

    out_name: str, optional, defaults to 'scatterplot'
        name of the output plot
    out_type: str, optional, defaults to 'pdf'
        type of the output plot
    out_directory: str, optional, defaults to ''
        directory path of the output plot
    
    Returns
    ----------
    None
    """

    df_data = pd.read_csv(df_filepath)
    df_filtered = df_data.loc[df_data[hue_column].isin(hue_order)]
    
    # make plot
    _, ax = plt.subplots(**subplots_kws)
    sns.scatterplot(data=df_filtered, 
                    ax=ax, 
                    x=condition1, y=condition2, 
                    hue=hue_column, hue_order=hue_order, palette=palette, 
                    **scatterplot_kws,
                    )
    
    # adjust x and y axis limits
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    # set labels
    plt.xlabel(xlab) # set x-axis label
    plt.ylabel(ylab) # set y-axis label
    plt.tight_layout()

    # save to pdf and close
    if savefig: 
        output_path = out_directory + condition1 + condition2 + '_' + out_name + '.' + out_type
        plt.savefig(output_path, format=out_type)
    plt.show()
    plt.close()

# python3 -m be_scan plot_corr_scatterplot -df '../../../Downloads/NZL10196_v9_comparisons.csv' -c1 'd3-neg' -c2 'd9-pos' -hue 'Mut_type'
