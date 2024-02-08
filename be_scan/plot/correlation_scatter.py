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
import scipy.stats as stats

from be_scan.plot._annotating_ import list_muttypes, color_list

def plot_corr_scatterplot(df_filepath, condition1, condition2, 

    filter_val=False, val_cols=[], val_min=None, # filter out unwanted quantitative params
    filter_params=False, params_cols=[], params_conditions=[], # filter out unwanted categorical params
    include_hue=False, hue_col='Mut_type', hue_order=list_muttypes, palette=color_list, # color params
    savefig=True, out_directory='', out_name='correlation_scatterplot', out_type='png', show=True,

    jointplot_kws={'alpha':0.8, 'linewidth':1, 
                    'edgecolor':'black', 's':25}, 
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
    include_hue: bool, optional, default to False
        whether or not to color points by a variable, 
        will also restrict points plotted to only the hue_order values listed
    hue_col: str, optional, defaults to 'Mut_type'
        the categorial dimension of the data, name of .csv data column
    hue_order: list of str, optional, defaults to a preset list of column names
        a list of categorial variables in hue_col
    palette: list of str, optional, defaults to a preset list of colors from ColorBrewer2
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
    show : bool, optional, defaults to True
        whether or not to show the plot

    jointplot_kws: dict, optional, defaults to 
        {'alpha':0.8, 'linewidth':1, 
        'edgecolor':'black', 's':25}
        input params for sns.jointplot() 
        https://seaborn.pydata.org/generated/seaborn.jointplot.html

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
    if filter_params: 
        assert isinstance(params_cols, list), "check param: params_cols"
        assert isinstance(params_conditions, list), "check param: params_conditions"
    if include_hue: 
        assert hue_col in df_data.columns.tolist(), "check param: hue_col"
        assert isinstance(hue_order, list) and len(hue_order) > 0, "check param: hue_order"
        assert isinstance(palette, list) and len(palette) > 0, "check param: palette"
    assert isinstance(jointplot_kws, dict), "check param: jointplot_kws"
    
    # apply 2 layers of filters
    if filter_params: 
        for col, conds in zip(params_cols, params_conditions): 
            df_data = df_data.loc[df_data[col].isin(conds)]
    if filter_val: 
        for col in val_cols: 
            df_data = df_data[df_data[col] > val_min]
            
    # set plot params, if there is hue add params for plotting
    baseline_params = {'data':df_data, 'x':condition1, 'y':condition2}
    if include_hue: 
        df_data = df_data.loc[df_data[hue_col].isin(hue_order)]
        baseline_params = {**baseline_params, 
                           **{'hue':hue_col, 'hue_order':hue_order, 'palette':palette}}
    
    # plot with params as input
    sns.jointplot(**baseline_params,
                  **jointplot_kws,
                 )
    plt.tight_layout()
    # set legend on the far right
    ax = plt.gca()
    ax.legend(loc='center left', bbox_to_anchor=(1.3, 0.5))

    # calculate some stats and output them
    slope, intercept, r_value, p_value, std_err = stats.linregress(df_data[condition1], 
                                                                   df_data[condition2])
    print("R: {0} (p-value {1})".format(r_value, p_value))
    print("R2: {0}".format(r_value**2))

    # save to pdf and close
    path = Path.cwd()
    if savefig: 
        outpath = path / out_directory
        out = condition1 + condition2 + '_' + out_name + '.' + out_type
        plt.savefig(outpath / out, format=out_type)
    if show: 
        plt.show()
    plt.close()
