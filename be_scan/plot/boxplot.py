"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function plots chosen guides by plot_column categories to show the distribution of guides}
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import matplotlib as mpl
import plotly.express as px
from plotly.subplots import make_subplots

from be_scan.plot._annotating_ import *
from be_scan.plot._shared_ import (
    annotate_points, apply_filters, control_stdev, draw_control_lines,
    load_dataframe, normalize_to_controls, resolve_hue, save_and_show,
    save_and_show_interactive, select_annotated, split_highlight, style_axes,
)

def boxplot(
    df_filepath, # DF Filename or DF #
    comparisons, # A List of Y-Axis Values #
    plot_column, # The Column of plot_conditions #
    plot_conditions, # A List of Boxes in the Plot #

    xlab='', ylab='Log2 Fold Change', # Scatterplot Labels #
    neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # Negative Ctrl Parameters #
    savefig=True, show=True, out_name='boxes', out_type='png', out_dir='', # Output Parameters #

    # Style Params #
    subplots_kws={},
    boxplot_kws = {'saturation':1, 'fliersize':4, 'width':0.4,
                   'flierprops':{'marker':'o', 'mec':'black', 'lw':1, 'alpha':0.8}},
    axhline_kws = {'color':'k', 'ls':'--', 'lw':1},

    # Added in 3.0.0 for parity with the other be_scan.plot functions. #
    # Appended so existing positional and keyword calls are unaffected, and #
    # every default leaves behavior exactly as it was. #
    title=None, dpi=300, n_sd=2,
    include_hue=False, hue_col='', pal='pastel', # Color Parameters #
    filter_val=False, val_cols=[], val_min=0.0, # Filter Out Values #
    filter_params=False, params_cols=[], params_conditions=[['Missense', 'Silent', 'Mixed', 'Nonsense']], # Filter Out Categories #
    annot=False, annot_label='', annot_top=None, annot_cutoff=None, annot_abs=None, # Annotate Outliers #
    highlight=False, highlight_col='', highlight_vals=[], # Highlight a Subset #
    xlim_kws=None, ylim_kws=None,
    label_text_kws={'ha':'right', 'va':'bottom', 'fontsize':6},
    stripplot_kws={'size':3, 'linewidth':0.5, 'edgecolor':'black', 'color':'red'},
    ):

    """[Summary]
    This function takes in a dataframe from count_reads, performs normalization,
    and then plots chosen (ex control) guides by plot_column categories,
    to show the distribution of categories of guides

    Parameters
    ------------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    comparisons : list of str, required
        list of comparisons that correspond to columns of data
    plot_column : str, required
        column of .csv, typically domain or mutation type
    plot_conditions : list of str, required
        category names of plot_column

    xlab : str, optional, defaults to ''
        name of x-axis label
    ylab : str, optional, defaults to 'Log2 Fold Change'
        name of y-axis label

    neg_ctrl : bool, optional, defaults to False
        whether or not to calulate negative control for normalization and line drawing
    neg_ctrl_col : str, optional, defaults to 'CtoT_win_overlap'
        column of .csv which correspond to normalization control
    neg_ctrl_conditions : list of str, optional, defaults to ['Intron', 'Control']
        name of categories of neg_ctrl_col to normalize dataframe

    savefig : boolean, optional, defaults to True
        option of saving figure to output or not
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
        path to output directory

    subplots_kws : dict, optional, defaults to
        {}
        input params for plt.subplots()
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
    boxplot_kws : dict, optional, defaults to
        {'saturation':1, 'fliersize':4, 'width':0.4,
        'flierprops':{'marker':'o', 'mec':'black', 'lw':1, 'alpha':0.8}}
        input params for sns.boxplot()
        https://seaborn.pydata.org/generated/seaborn.boxplot.html
    axhline_kws : dict, optional, defaults to
        {'color':'k', 'ls':'--', 'lw':1}
        input params for plt.axhline()
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axhline.html

    Returns
    ------------
    fig, axes : the matplotlib figure and its list of axes
    """
    # LOAD DATAFRAME #
    df_data = load_dataframe(df_filepath)

    # Normalize Data to Intergenic Ctrls #
    df_data, list_ctrlstats = normalize_to_controls(
        df_data, comparisons, neg_ctrl, neg_ctrl_col, neg_ctrl_conditions)

    # Filter Data Before Plotting #
    df_data = apply_filters(
        df_data, filter_val, val_cols, val_min,
        filter_params, params_cols, params_conditions)

    df_data = df_data.loc[df_data[plot_column].isin(plot_conditions)].copy()

    # Add Hue #
    hue_params = resolve_hue(df_data, include_hue, hue_col, pal)

    # Autoannotate the Hits #
    if annot:
        assert isinstance(annot_label, str), "check param: annot_label"
        assert annot_label in df_data.columns.tolist(), "check param: annot_label"

    # Setup Subplots #
    fig, axes = plt.subplots(
        nrows=len(comparisons), ncols=1,
        **subplots_kws)
    if len(comparisons) == 1: axes = [axes]

    for ax, comp in zip(axes, comparisons):
        # Plot for Every Comparison #
        sns.boxplot(data=df_data, ax=ax, x=plot_column, y=comp,
                    **hue_params, **boxplot_kws)
        plt.setp(ax.artists, edgecolor='black')
        plt.setp(ax.lines, color='black')

        # Overlay a Highlighted Subset #
        hits, _ = split_highlight(df_data, highlight, highlight_col, highlight_vals)
        if hits is not None and not hits.empty:
            sns.stripplot(data=hits, ax=ax, x=plot_column, y=comp, **stripplot_kws)

        # Overlay +/- n sd as Black Dashed Lines #
        draw_control_lines(ax, control_stdev(list_ctrlstats, comp),
                           n_sd=n_sd, **axhline_kws)

        # Annotate Hits #
        if annot:
            rows = select_annotated(df_data, comp, annot_top, annot_cutoff, annot_abs)
            annotate_points(ax, rows, plot_column, comp, annot_label,
                            **label_text_kws)

        # Set Title and Axes #
        # ylim_kws overrides the data-derived range when supplied
        if ylim_kws is None:
            ax.set_ylim(np.floor(df_data[comp].min()), np.ceil(df_data[comp].max()))
        style_axes(ax, xlab=xlab, ylab=ylab,
                   title=comp if title is None else title,
                   xlim_kws=xlim_kws, ylim_kws=ylim_kws, hide_spines=())

        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_horizontalalignment('right')

    # Save File #
    save_and_show(fig, savefig=savefig, show=show, out_dir=out_dir,
                  out_name=out_name, out_type=out_type, dpi=dpi)
    return fig, axes

def interactive_boxplot(
    df_filepath, # DF Filename or DF #
    comparisons, # A List of Y-Axis Values #
    plot_column, # The Column of plot_conditions #
    plot_conditions, # A List of Boxes in the Plot #

    figsize=(500, 300),
    xlab='', ylab='Log2 Fold Change', # Scatterplot Labels #
    neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # Negative Ctrl Parameters #
    savefig=True, show=True, out_name='boxes', out_type='png', out_dir='', # Output Parameters #

    # Added in 3.0.0 for parity with the other interactive plots. #
    title='Box Plots - Comparsions',
    include_hue=False, hue_col='', pal='pastel', # Color Parameters #
    filter_val=False, val_cols=[], val_min=0.0, # Filter Out Values #
    filter_params=False, params_cols=[], params_conditions=[['Missense', 'Silent', 'Mixed', 'Nonsense']], # Filter Out Categories #
    ):

    """[Summary]
    This function takes in a dataframe from count_reads, performs normalization,
    and then plots chosen (ex control) guides by plot_column categories,
    to show the distribution of categories of guides

    Parameters
    ------------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    comparisons : list of str, required
        list of comparisons that correspond to columns of data
    plot_column : str, required
        column of .csv, typically domain or mutation type
    plot_conditions : list of str, required
        category names of plot_column

    figsize : tuple, optional, defaults to (500, 300)
        (width, height) of figure in pixels
    xlab : str, optional, defaults to ''
        name of x-axis label
    ylab : str, optional, defaults to 'Log2 Fold Change'
        name of y-axis label

    neg_ctrl : bool, optional, defaults to False
        whether or not to calulate negative control for normalization and line drawing
    neg_ctrl_col : str, optional, defaults to 'CtoT_win_overlap'
        column of .csv which correspond to normalization control
    neg_ctrl_conditions : list of str, optional, defaults to ['Intron', 'Control']
        name of categories of neg_ctrl_col to normalize dataframe

    savefig : boolean, optional, defaults to True
        option of saving figure to output or not
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_dir : str, optional, defaults to ''
        path to output directory

    Returns
    ------------
    """
    # LOAD DATAFRAME #
    df_data = load_dataframe(df_filepath)

    # Normalize Data to Intergenic Ctrls #
    df_data, list_ctrlstats = normalize_to_controls(
        df_data, comparisons, neg_ctrl, neg_ctrl_col, neg_ctrl_conditions)

    # Filter Data Before Plotting #
    df_data = apply_filters(
        df_data, filter_val, val_cols, val_min,
        filter_params, params_cols, params_conditions)

    df_data = df_data.loc[df_data[plot_column].isin(plot_conditions)].copy()

    # Setup Subplots #
    fig = make_subplots(
        rows=len(comparisons), cols=1,
        subplot_titles=comparisons,
        )

    for i, comp in enumerate(comparisons):
        box_fig = px.box(df_data, x=plot_column, y=comp,
                         color=hue_col if include_hue else None)
        for trace in box_fig.data:
            fig.add_trace(trace, row=i+1, col=1)

        # Overlay neg ctrl avg +/- 2 SD as black dashed lines
        stdev = control_stdev(list_ctrlstats, comp)
        if stdev is not None:
            fig.add_hline(y=2*stdev, line=dict(dash="dash", color="black"), row=i+1, col=1)
            fig.add_hline(y=-2*stdev, line=dict(dash="dash", color="black"), row=i+1, col=1)

        # Adjust Y Limits #
        fig.update_yaxes(
            range=[np.floor(df_data[comp].min()), np.ceil(df_data[comp].max())],
            title_text=ylab, row=i+1, col=1
        )

    # Update Layout #
    fig.update_layout(
        title=title,
        xaxis_title=xlab, yaxis_title=ylab,
        width=figsize[0], height=figsize[1] * len(comparisons),  # Adjust height dynamically
        showlegend=include_hue,
    )

    # Save File #
    save_and_show_interactive(fig, savefig=savefig, show=show,
                              out_dir=out_dir, out_name=out_name)
    return fig
