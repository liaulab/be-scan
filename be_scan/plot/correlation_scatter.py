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

def plot_corr_scatterplot(df_filepath, 
                          condition1, condition2, 

                          filter=False, filter_col=[], filter_conditions=[], # filter out unwanted data params
                          hue=False, hue_column='Mut_type', hue_order=list_muttypes, palette=color_list, # color params
                          savefig=True, out_directory='', out_name='correlation_scatterplot', out_type='pdf', 

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

    jointplot_kws: dict, optional, defaults to {'alpha':0.8, 'linewidth':1, 'edgecolor':'black', 's':25}
        input params for sns.jointplot()

    Returns
    ----------
    None
    """

    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)
    
    if filter: 
        for col, conds in zip(filter_col, filter_conditions): 
            df_data = df_data.loc[df_data[col].isin(conds)]
            
    # make plot
    baseline_params = {'data':df_data, 'x':condition1, 'y':condition2}
    if hue: 
        df_data = df_data.loc[df_data[hue_column].isin(hue_order)]
        baseline_params = {**baseline_params, 
                           **{'hue':hue_column, 'hue_order':hue_order, 'palette':palette}}
        
    sns.jointplot(**baseline_params,
                  **jointplot_kws,
                 )
    plt.tight_layout()
    ax = plt.gca()
    ax.legend(loc='center left', bbox_to_anchor=(1.3, 0.5))

    slope, intercept, r_value, p_value, std_err = stats.linregress(df_data[condition1], df_data[condition2])
    print("R: {0} (p-value {1})".format(r_value, p_value))
    print("R2: {0}".format(r_value**2))

    # save to pdf and close
    path = Path.cwd()
    if savefig: 
        outpath = path / out_directory
        out = condition1 + condition2 + '_' + out_name + '.' + out_type
        plt.savefig(outpath / out, format=out_type)
    plt.show()
    plt.close()
