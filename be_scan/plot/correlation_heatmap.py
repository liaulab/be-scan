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
import plotly.graph_objects as go

def corr_heatmap(
    df_filepath, # DF Filename or DF #
    comparisons, # A List of Y-Axis Values #

    corr_type='spearman',
    xlab='', ylab='', title='Spearman Correlation Heatmap', # Labels #
    savefig=True, show=True, out_dir='', out_name='correlation_heatmap', out_type='png', # Output Parameters #

    # Style Params #
    heatmap_kws={'center':0, 'linewidth':0.5, 'cmap':'coolwarm',
                 'square':True, 'cbar_kws':{"shrink": 0.5}, 'annot':True},
    subplots_kws = {},
    ):

    """[Summary]
    This function takes in a dataframe from count_reads, and plots
    a scatterplot showing correlation between two given conditions

    Parameters
    ------------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    comparisons : list of str, required
        list of comparisons that correspond to columns of .csv data

    corr_type : str, optional, defaults to 'spearman'
        type of correlation, refer to https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.corr.html
    xlab : str, optional, defaults to ''
        name of the x-axis label
    ylab : str, optional, defaults to ''
        name of the y-axis label
    title : str, optional, defaults to 'Spearman Correlation Heatmap'
        name of title label

    savefig : bool, optional, defaults to True
        whether or not to save the figure
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_dir : str, optional, defaults to ''
        path to output directory
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output

    heatmap_kws : dict, optional, defaults to
        {'center':0, 'linewidth':0.5, 'cmap':'coolwarm',
        'square':True, 'cbar_kws':{"shrink": 0.5}, 'annot':True}
        input params for sns.heatmap()
        https://seaborn.pydata.org/generated/seaborn.heatmap.html
    subplots_kws : dict, optional, defaults to
        {}
        input params for plt.subplots()
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html

    Returns
    ------------
    """
    # LOAD DATAFRAME #
    if isinstance(df_filepath, pd.DataFrame):
        df_data = df_filepath.copy()
    elif isinstance(df_filepath, (str, Path)):
        df_filepath = Path(df_filepath)
        df_data = pd.read_csv(df_filepath)
    else:
        raise ValueError("df_filepath must be a pandas DataFrame or a path string/Path.")

    sns.set_style('ticks')

    # Compute Correlation Matrix #
    df_comp = df_data[comparisons].copy()
    df_corr = df_comp.corr(method=corr_type)

    # Setup Subplots #
    _, ax = plt.subplots(**subplots_kws)
    ax = sns.heatmap(df_corr, **heatmap_kws)

    # Adjustments and Labels #
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    plt.title(title) ; plt.ylabel(xlab) ; plt.xlabel(ylab)
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.yticks(rotation=0, horizontalalignment='right')
    plt.tight_layout()

    # Save File #
    outpath = Path(out_dir)
    if savefig:
        out = f'{out_name}.{out_type}'
        plt.savefig(outpath / out, format=out_type, dpi=300)
    if show: plt.show()
    plt.close()

def interactive_corr_heatmap(
    df_filepath, # DF Filename or DF #
    comparisons, # A List of Y-Axis Values #

    corr_type='spearman',
    figsize=(500,500),
    xlab='', ylab='', title='Spearman Correlation Heatmap', # Labels #
    savefig=True, show=True, out_dir='', out_name='correlation_heatmap', # Output Parameters #
    ):

    """[Summary]
    This function takes in a dataframe from count_reads, and plots
    a scatterplot showing correlation between two given conditions

    Parameters
    ------------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    comparisons : list of str, required
        list of comparisons that correspond to columns of .csv data

    corr_type : str, optional, defaults to 'spearman'
        type of correlation, refer to https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.corr.html
    figsize : tuple, optional, defaults to (500, 300)
        (width, height) of figure in pixels
    xlab : str, optional, defaults to ''
        name of the x-axis label
    ylab : str, optional, defaults to ''
        name of the y-axis label
    title : str, optional, defaults to 'Spearman Correlation Heatmap'
        name of title label

    savefig : bool, optional, defaults to True
        whether or not to save the figure
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_dir : str, optional, defaults to ''
        path to output directory
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output

    Returns
    ------------
    """
    # LOAD DATAFRAME #
    if isinstance(df_filepath, pd.DataFrame):
        df_data = df_filepath.copy()
    elif isinstance(df_filepath, (str, Path)):
        df_filepath = Path(df_filepath)
        df_data = pd.read_csv(df_filepath)
    else:
        raise ValueError("df_filepath must be a pandas DataFrame or a path string/Path.")

    sns.set_style('ticks')

    # Compute Correlation Matrix #
    df_comp = df_data[comparisons].copy()
    df_corr = df_comp.corr(method=corr_type)

    fig = go.Figure()
    # Add Heatmap #
    fig.add_trace(go.Heatmap(
        z=df_corr.values, x=df_corr.columns, y=df_corr.index,
        colorscale="rdbu_r",
        colorbar=dict(title="Correlation"),
        zmin=-1, zmax=1  # Ensure scale is between -1 and 1 for correlation
    ))
    # Adjustments and Labels #
    fig.update_layout(
        title=title,
        width=figsize[0], height=figsize[1],
        xaxis=dict(title=xlab, tickangle=-45), yaxis=dict(title=ylab, tickangle=0),
        margin=dict(l=50, r=50, t=50, b=50),
    )

    # Save File #
    outpath = Path(out_dir)
    if savefig:
        fig.write_html(outpath / f"{out_name}.html")
    if show: fig.show()
