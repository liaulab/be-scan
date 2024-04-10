"""
Author: Calvin XiaoYang Hu
Date: 240314

{Description: This function highlights certain guides on a scale of their logFC values, 
              and shows the distribution of logFC values}
"""

import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def sgrna_enrichment(df_filepath, comparisons, 
                     
     highlight=False, highlight_col='sgRNA_ID', highlight_vals=['sgRNA_0'], 
     density=True, 
     
     savefig=True, out_name='scatterplot', out_type='png', out_directory='', show=True, # output params
                     
     rugplot_kws={'height':0.25, 'color':'black', 'linewidth':2, 'alpha':0.05}, 
     kdeplot_kws={'color':'black', 'linewidth':1, 'alpha':0.5}, 
     highlight_rugplot_kws={'height':0.25, 'color':'red', 'linewidth':5, 'alpha':0.6}
    ):
    
    """[Summary]
    This function takes in a dataframe from analysis and then 
    plots the data for each condition to highlight specifically enriched guides

    Parameters
    ----------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    comparisons : list of str, required
        list of comparisons that correspond to columns of data
        each comparison is a plot
    
    highlight : bool, optional, defaults to False
        whether or not to highlight specific guides
    highlight_col : str, optional, defaults to ''
        name of a column to filter for guides
    highlight_conditions : list of str, optional, defaults to []
        names of values in highlight_col to highlight
    density: bool, optional, defaults to True
        whether or not to include density portion of the rugplot
        
    savefig: bool, optional, defaults to True
        whether or not to save the figure
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
        path to output directory
    show : bool, optional, defaults to True
        whether or not to show the plot
        
    rugplot_kws: dict, optional, defaults to
        {'height':0.25, 'color':'black', 'linewidth':2, 'alpha':0.05}
        additional input of seaborn.rugplot()
        https://seaborn.pydata.org/generated/seaborn.rugplot.html
    kdeplot_kws: dict, optional, defaults to
        {'color':'black', 'linewidth':1, 'alpha':0.5}
        additional input of seaborn.kdeplot_kws()
        https://seaborn.pydata.org/generated/seaborn.kdeplot.html
    highlight_rugplot_kws: dict, optional, defaults to
        {'height':0.25, 'color':'red', 'linewidth':5, 'alpha':0.6}
        additional input of seaborn.rugplot() for highlighted guides
        https://seaborn.pydata.org/generated/seaborn.rugplot.html

    Returns
    ----------
    None
    """
    
    # style
    mpl.rcParams.update({'font.size': 10})
    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)
    
    # set up same amount of plots as comparisons
    # figsize is dependent on number of plots
    fig, ax = plt.subplots(nrows=len(comparisons), ncols=1, 
                           figsize=(30, 5*len(comparisons)))

    if highlight: 
        assert isinstance(highlight_col, str), 'check param: highlight_col'
        assert isinstance(highlight_vals, list), 'check param: highlight_vals'
    assert isinstance(rugplot_kws, dict), "check param: rugplot_kws"
    assert isinstance(kdeplot_kws, dict), "check param: kdeplot_kws"
    assert isinstance(highlight_rugplot_kws, dict), "check param: highlight_rugplot_kws"
    
    # add a density plot
    if len(comparisons) == 1: 
        ax = [ax]
    if density: 
        for i, val in enumerate(comparisons): 
            sns.kdeplot(data=df_data, x=val, 
                        ax=ax[i], **kdeplot_kws)
    # add a rugplot to show density of guides across enrichment values
    for i, val in enumerate(comparisons): 
        sns.rugplot(data=df_data, x=val, 
                    ax=ax[i], **rugplot_kws)
        ax[i].set_title(f"sgRNA Enrichment for {val}")
        # add option to highlight certain guides
        if highlight: 
            highlight_df = df_data[df_data[highlight_col].isin(highlight_vals)]
            sns.rugplot(data=highlight_df, x=val, 
                        ax=ax[i], **highlight_rugplot_kws)

    plt.tight_layout()
    ### title and labels
    
    # Save to pdf
    path = Path.cwd()
    if savefig:
        outpath = path / out_directory
        out = out_name + '.' + out_type
        plt.savefig(outpath / out, format=out_type)
    if show: 
        plt.show()
    plt.close()
    