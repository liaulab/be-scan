"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function plots a heatmap showing correlation between all conditions}
"""

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def corr_heatmap(df_filepath, comparisons, 

    corr_type='spearman', 
    xlab='', ylab='', title='Spearman Correlation Heatmap', # figure related params
    savefig=True, show=True, out_dir='', out_name='correlation_heatmap', out_type='png', # output related params

    # style params
    heatmap_kws={'center':0, 'linewidth':0.5, 'cmap':'coolwarm', 
                 'square':True, 'cbar_kws':{"shrink": 0.5}, 'annot':True},
    subplots_kws = {'figsize':(4,4)}, 
    ): 
    
    """[Summary]
    This function takes in a dataframe from count_reads, and plots
    a scatterplot showing correlation between two given conditions
    
    Parameters
    ----------
    df_filepath: str, required
        filepath to .csv data generated from count_reads
    comparisons: list of str, required
        list of comparisons that correspond to columns of .csv data

    corr_type : str, optional, defaults to 'spearman'
        type of correlation, refer to https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.corr.html
    xlab : str, optional, defaults to ''
        name of the x-axis label
    ylab : str, optional, defaults to ''
        name of the y-axis label
    title : str, optional, defaults to 'Spearman Correlation Heatmap'
        name of title label
    savefig: bool, optional, defaults to True
        whether or not to save the figure
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_dir : str, optional, defaults to ''
        path to output directory
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    
    heatmap_kws: dict, optional, defaults to 
        {'center':0, 'linewidth':0.5, 'cmap':'coolwarm', 
        'square':True, 'cbar_kws':{"shrink": 0.5}, 'annot':True}
        input params for sns.heatmap() 
        https://seaborn.pydata.org/generated/seaborn.heatmap.html
    subplots_kws: dict, optional, defaults to 
        {'figsize':(4,4)}
        input params for plt.subplots() 
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html

    Returns
    ----------
    None
    """

    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)

    # style
    mpl.rcParams.update({'font.size': 10})
    sns.set_style('ticks')

    # Compute correlation matrix
    df_comp = df_data[comparisons].copy()
    df_corr = df_comp.corr(method=corr_type)

    # Set up the matplotlib figure
    _, ax = plt.subplots(**subplots_kws)
    ax = sns.heatmap(df_corr, **heatmap_kws)
    
    # show frame
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    # adjustments and labels
    plt.title(title) ; plt.ylabel(xlab) ; plt.xlabel(ylab)
    # rotate axis labels
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.yticks(rotation=0, horizontalalignment='right')
    plt.tight_layout()

    # save pdf and close everything
    outpath = Path(out_dir)
    if savefig: 
        out_name = f'{out_name}.{out_type}'
        plt.savefig(outpath / out_name, format=out_type, dpi=300)
    if show: plt.show()
    plt.close()
