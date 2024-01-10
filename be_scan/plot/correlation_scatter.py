"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function plots a scatterplot showing correlation between two given conditions}
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from be_scan.plot._annotating_ import list_muttypes, color_list

def plot_corr_scatterplot(df_filepath, 
                          condition1, condition2, 

                          hue=False, hue_column='Mut_type', hue_order=list_muttypes, palette=color_list, # color params
                          xlab='cond1 score', ylab='cond2 score', 
                          savefig=True, out_directory='', out_name='correlation_scatterplot', out_type='pdf', 

                          xlim_kws={'xmin':None, 'xmax':None}, ylim_kws={'ymin':None, 'ymax':None},
                          scatterplot_kws={'alpha':0.8, 'linewidth':1, 
                                           'edgecolor':'black', 's':25}, 
                          subplots_kws = {'figsize':(4.5, 4)},
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

    hue: bool, optional, default to False
        whether or not to color points by a variable, will also restrict points plotted to only the hue_order values listed
    hue_column: str, optional, defaults to 'Mut_type'
        the categorial dimension of the data, name of .csv data column
    hue_order: list of str, optional, defaults to list_muttypes in _annotating_.py a preset list of column names
        a list of categorial variables in hue_column
    palette: list of str, optional, defaults to color_list in _annotating_.py a preset list of colors from ColorBrewer2
        a list of colors which correspond to hue_order
    xlab: str, optional, defaults to 'cond1 score'
        x-axis label
    ylab: str, optional, defaults to 'cond2 score'
        y-axis label
        
    savefig: bool, optional, defaults to True
        whether or not to save the figure
    out_name: str, optional, defaults to 'scatterplot'
        name of the output plot
    out_type: str, optional, defaults to 'pdf'
        type of the output plot
    out_directory: str, optional, defaults to ''
        directory path of the output plot

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

    df_data = pd.read_csv(df_filepath)
    
    # make plot
    _, ax = plt.subplots(**subplots_kws)
    baseline_params = {'data':df_data, 'ax':ax, 'x':condition1, 'y':condition2}
    if hue: 
        df_data = df_data.loc[df_data[hue_column].isin(hue_order)]
        sns.scatterplot(**baseline_params,
                        hue=hue_column, hue_order=hue_order, palette=palette, 
                        **scatterplot_kws,
                        )
    else: 
        sns.scatterplot(**baseline_params,
                        **scatterplot_kws,
                        )
    
    # adjust x and y axis limits
    ax.set_xlim(**xlim_kws)
    ax.set_ylim(**ylim_kws)
    # set labels
    plt.xlabel(xlab) # set x-axis label
    plt.ylabel(ylab) # set y-axis label
    plt.tight_layout()

    # save to pdf and close
    if savefig: 
        out = condition1 + condition2 + '_' + out_name + '.' + out_type
        output_path = Path(out_directory)
        plt.savefig(output_path / out, format=out_type)
    plt.show()
    plt.close()
