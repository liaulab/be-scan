"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function plots a heatmap showing correlation between all conditions}
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def plot_corr_heatmap(df_filepath, comparisons, 

    corr_type='spearman', 
    filter_val=False, val_cols=[], val_min=None, # filter out unwanted quantitative params
    filter_params=False, params_cols=[], params_conditions=[], # filter out unwanted categorical params
    xlab='', ylab='', title='Spearman Correlation Heatmap', # figure related params
    savefig=True, out_directory='', out_name='correlation_heatmap', out_type='png', show=True, # output related params

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
    xlab : str, optional, defaults to ''
        name of the x-axis label
    ylab : str, optional, defaults to ''
        name of the y-axis label
    title : str, optional, defaults to 'Spearman Correlation Heatmap'
        name of title label

    savefig: bool, optional, defaults to True
        whether or not to save the figure
    out_directory : str, optional, defaults to ''
        path to output directory
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    show : bool, optional, defaults to True
        whether or not to show the plot
    
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

    # check conflicting params and output for user
    if filter_val: 
        assert isinstance(val_min, float), "check param: val_min"
        assert isinstance(val_cols, list) and len(val_cols) > 0, "check param: val_cols"
        for vc in val_cols: 
            assert vc in df_data.columns.tolist(), "check param: val_cols"
    if filter_params: 
        assert isinstance(params_cols, list), "check param: params_cols"
        for pc in params_cols: 
            assert pc in df_data.columns.tolist(), "check param: val_cols"
        assert isinstance(params_conditions, list), "check param: params_conditions"
    assert isinstance(heatmap_kws, dict), "check param: heatmap_kws"
    assert isinstance(subplots_kws, dict), "check param: subplots_kws"
    assert corr_type in ["pearson", "kendall", "spearman"], "check param: corr_type"
    
    # apply 2 layers of filters
    if filter_params: 
        for col, conds in zip(params_cols, params_conditions): 
            df_data = df_data.loc[df_data[col].isin(conds)]
    if filter_val: 
        for col in val_cols: 
            df_data = df_data[df_data[col] > val_min]

    # Plotting parameters and variables
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
    plt.title(title)
    plt.ylabel(xlab)
    plt.xlabel(ylab)
    # rotate axis labels
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.yticks(rotation=0, horizontalalignment='right')
    plt.tight_layout()

    # save pdf and close everything
    path = Path.cwd()
    if savefig: 
        outpath = path / out_directory
        out_name = out_name + '.' + out_type
        plt.savefig(outpath / out_name, format=out_type)
    if show: 
        plt.show()
    plt.close()
