"""
Author: Calvin XiaoYang Hu
Date: 240314

{Description: This function highlights certain guides on a scale of their logFC values, 
              and shows the distribution of logFC values}
"""

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

from be_scan.plot._shared_ import (
    apply_filters, load_dataframe, normalize_to_controls, resolve_hue,
    save_and_show, split_highlight, style_axes,
)

def sgrna_enrichment(
    df_filepath, # DF Filename or DF #
    comparisons, # A List of Y-Axis Values #

    highlight=False, highlight_col='sgRNA_ID', highlight_vals=['sgRNA_0'],
    density=True,

    savefig=True, show=True, out_name='scatterplot', out_type='png', out_dir='', # Output Params #

    subplots_kws={},
    rugplot_kws={'height':0.25, 'color':'black', 'linewidth':2, 'alpha':0.05},
    kdeplot_kws={'color':'black', 'linewidth':1, 'alpha':0.5},
    highlight_rugplot_kws={'height':0.25, 'color':'red', 'linewidth':3, 'alpha':0.6},

    # Added in 3.0.0 for parity with the other be_scan.plot functions. #
    # Appended so existing calls are unaffected, and every default leaves #
    # behavior exactly as it was. #
    #
    # Note: out_name above defaults to 'scatterplot' rather than
    # 'sgrna_enrichment'. That is a copy-paste artifact, but correcting it
    # would silently change the output filename for existing callers, so it
    # is left alone deliberately. Pass out_name explicitly.
    xlab=None, ylab=None, title=None, dpi=300,
    neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # Negative Ctrl Parameters #
    filter_val=False, val_cols=[], val_min=0.0, # Filter Out Values #
    filter_params=False, params_cols=[], params_conditions=[['Missense', 'Silent', 'Mixed', 'Nonsense']], # Filter Out Categories #
    include_hue=False, hue_col='', pal='pastel', # Color Parameters #
    ):

    """[Summary]
    This function takes in a dataframe from analysis and then
    plots the data for each condition to highlight specifically enriched guides

    Parameters
    ------------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    comparisons : list of str, required
        list of comparisons that correspond to columns of .csv data

    highlight : bool, optional, defaults to False
        whether or not to highlight specific guides
    highlight_col : str, optional, defaults to ''
        name of a column to filter for guides
    highlight_conditions : list of str, optional, defaults to []
        names of values in highlight_col to highlight
    density: bool, optional, defaults to True
        whether or not to include density portion of the rugplot

    savefig : bool, optional, defaults to True
        whether or not to save the figure
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_dir : str, optional, defaults to ''
        path to output directory

    subplots_kws : dict, optional, defaults to {}
        input params for plt.subplots()
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
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
    ------------
    """
    # LOAD DATAFRAME #
    df_data = load_dataframe(df_filepath)

    # Normalize Data to Intergenic Ctrls #
    df_data, _ = normalize_to_controls(
        df_data, comparisons, neg_ctrl, neg_ctrl_col, neg_ctrl_conditions)

    # Filter Data Before Plotting #
    df_data = apply_filters(
        df_data, filter_val, val_cols, val_min,
        filter_params, params_cols, params_conditions)

    # Add Hue #
    hue_params = resolve_hue(df_data, include_hue, hue_col, pal)

    # Set Up Subplots #
    fig, ax = plt.subplots(
        nrows=len(comparisons), ncols=1,
        **subplots_kws)

    if highlight:
        assert isinstance(highlight_col, str) and highlight_col in df_data.columns, 'check param: highlight_col'
        assert isinstance(highlight_vals, list), 'check param: highlight_vals'

    # Add Density Plot #
    if len(comparisons) == 1:
        ax = [ax]
    if density:
        for i, val in enumerate(comparisons):
            sns.kdeplot(
                data=df_data, x=val, ax=ax[i],
                **hue_params, **kdeplot_kws)
    # Add Rugplot #
    for i, val in enumerate(comparisons):
        sns.rugplot(
            data=df_data, x=val, ax=ax[i],
            **rugplot_kws)
        # Highlight Particular Guides #
        if highlight:
            highlight_df, _ = split_highlight(
                df_data, highlight, highlight_col, highlight_vals)
            sns.rugplot(
                data=highlight_df, x=val, ax=ax[i],
                **highlight_rugplot_kws)
        style_axes(ax[i], xlab=xlab, ylab=ylab,
                   title=f"sgRNA Enrichment for {val}" if title is None else title,
                   hide_spines=())

    # Save File #
    save_and_show(fig, savefig=savefig, show=show, out_dir=out_dir,
                  out_name=out_name, out_type=out_type, dpi=dpi)
    return fig, ax
