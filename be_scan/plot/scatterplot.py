"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function performs normalization and then plots the data for each condition to reveal enriched guides}
"""

import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import warnings
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matplotlib.colors as mcolors

from be_scan.plot._annotating_ import *
from be_scan.plot._shared_ import (
    annotate_points, apply_filters, apply_xwindow, control_stdev,
    draw_control_lines, draw_domains, hue_color_map, load_dataframe,
    normalize_to_controls, resolve_hue, save_and_show,
    save_and_show_interactive, select_annotated, split_highlight, style_axes,
)


def _display(df):
    """
    Show the annotated rows, in or out of a notebook.

    This was a bare call to `display`, which only exists as a builtin inside
    IPython. Any script using annot_top or annot_abs raised NameError right
    after drawing the labels.
    """
    try:
        from IPython.display import display as _ipy_display
    except ImportError:
        print(df.to_string())
        return
    try:
        from IPython import get_ipython
        if get_ipython() is None:
            print(df.to_string())
            return
    except ImportError:
        print(df.to_string())
        return
    _ipy_display(df)


def scatterplot(
    df_filepath, # DF Filename or DF #
    comparisons, # A List of Y-Axis Values #
    x_column, # X-Axis #

    xlab='Amino Acid Position', ylab='sgRNA Score', # Scatterplot Labels #
    xwindow=None, # Filter Data on xwindow #
    include_hue=False, hue_col='', pal='pastel', # Color Parameters #
    neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # Negative Ctrl Parameters #
    domains=[], domains_alpha=0.25, domains_color='lightblue', # Draw Domains #
    savefig=True, show=True, out_name='scatterplot', out_type='png', out_dir='', # Output Parameters #

    filter_val=False, val_cols=[], val_min=0.0, # Filter Out Values #
    filter_params=False, params_cols=[], params_conditions=[['Missense', 'Silent', 'Mixed', 'Nonsense']], # Filter Out Categories #
    annot=False, annot_label='', annot_top=None, annot_cutoff=None, annot_abs=None, # Annotate Outliers #

    # Style Params #
    subplots_kws={}, 
    xlim_kws={'xmin':None, 'xmax':None}, ylim_kws={'ymin':None, 'ymax':None},
    scatterplot_kws={'alpha': 0.85, 'linewidth': 0.8, 's': 30},
    axhline_kws={'color':'k', 'ls':'--', 'lw':1},
    label_text_kws={'ha':'right', 'va':'bottom', 'fontsize':6},

    # Added in 3.0.0 for parity with the other be_scan.plot functions. #
    # Appended so existing calls are unaffected, and every default leaves #
    # behavior exactly as it was. #
    title=None, dpi=300, n_sd=2, annot_offset=0.2,
    highlight=False, highlight_col='', highlight_vals=[], # Highlight a Subset #
    highlight_kws={'s':60, 'linewidth':1, 'edgecolor':'black', 'color':'red', 'zorder':10},
    ):

    """[Summary]
    This function takes in a dataframe from analysis and then
    plots the data for each condition to reveal which guides are enriched

    Parameters
    ------------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    x_column : str, required
        column of .csv, typically amino acid position
    comparisons : list of str, required
        list of comparisons that correspond to columns of data

    xlab : str, optional, defaults to 'Amino Acid Position'
        name of x-axis label
    ylab : str, optional, defaults to 'sgRNA Score'
        name of y-axis label
    xwindow : list, optional, defaults to []
        an inclusive range to restrict the x axis value to display

    include_hue: bool, optional, default to False
        whether or not to color points by a variable,
        will also restrict points plotted to only the hue_order values listed
    hue_col: str, optional, defaults to 'CtoT_muttype'
        the categorial dimension of the data, name of .csv data column
    palette: list of str, optional, defaults to a preset list of colors from ColorBrewer2
        a list of colors which correspond to hue_order

    filter_val : bool, optional, defaults to False
        whether or not to exclude a subset of data from plotting by a minimum value
        default purpose is to filter out all intron and exon/intron guides
    val_cols : list of str, optional,
        defaults to ['CtoT_muttypes']
        names of columns to filter dataframe for plotting
    val_min : int, optional, defaults to 0.0
        the minimum value by which to filter rows by val_cols

    filter_params : bool, optional, defaults to False
        whether or not to exclude a subset of data from plotting by categorical params
        default purpose is to filter to keep only Missense, Nonsense, Silent, Mixed guides
    params_cols : list of str, optional,
        defaults to ['CtoT_muttype']
        names of column to filter dataframe for plotting
    params_conditions : list of lists of str, optional,
        defaults to [['Missense', 'Silent', 'Mixed', 'Nonsense']]
        names of categories of filter_col to filter dataframe

    annot : bool, optional, defaults to False
        Whether or not to annot points
    annot_label : string, optional, defaults to 'CtoT_mutations'
        The column of the label in the dataframe
    annot_top : int, optional, defaults to 10
        The top n scoring points will be labeled
    annot_cutoff : float, optional, defaults to None
        The absolute value cutoff for which points will be labeled
    annot_abs : int, optional, defaults to None
        The absolute value cutoff for which points will be labeled

    domains: list, optional, defaults to []
        a list of domains with each item as {'start':n, 'end':m}
    domains_alpha: float, optional, defaults to 0.25
        alpha value of domains
    domains_color: str, optional, defaults to 'lightblue'
        color of domains

    neg_ctrl : bool, optional, defaults to False
        whether or not to calulate negative control for normalization and line drawing
    neg_ctrl_col : str, optional, defaults to 'CtoT_win_overlap'
        column of .csv which correspond to normalization control
    neg_ctrl_conditions : list of str, optional, defaults to ['Intron', 'Control']
        name of categories of neg_ctrl_col to normalize dataframe

    savefig : bool, optional, defaults to True
        whether or not to save the figure
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'png'
        file type of figure output
    out_dir : str, optional, defaults to ''
        path to output directory

    xlim_kws : dict, optional, defaults to
        {'xmin':None, 'xmax':None}
        input params for ax.set_xlim()
        https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_xlim.html
    ylim_kws : dict, optional, defaults to
        {'ymin':None, 'ymax':None}
        input params for ax.set_ylim()
        https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.set_ylim.html
    scatterplot_kws : dict, optional, defaults to
        {'alpha':0.8, 'linewidth':1.0, 'edgecolor':'black', 's':25}
        input params for sns.scatterplot()
        https://seaborn.pydata.org/generated/seaborn.scatterplot.html
    subplots_kws : dict, optional, defaults to 
        {}
        input params for plt.subplots()
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
    axhline_kws : dict, optional, defaults to
        {'color':'k', 'ls':'--', 'lw':1}
        input params for plt.axhline()
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axhline.html

    Returns
    ------------
    """

    # LOAD DATAFRAME #
    df_data = load_dataframe(df_filepath)

    # Normalize Data to Intergenic Ctrls #
    df_data, list_ctrlstats = normalize_to_controls(
        df_data, comparisons, neg_ctrl, neg_ctrl_col, neg_ctrl_conditions)

    # Filter Data Before Plotting #
    df_data = apply_xwindow(df_data, x_column, xwindow)

    # Filter by Filter Params #
    # applied before the hue mapping is built, so the legend only lists
    # categories that survive the filters
    df_data = apply_filters(
        df_data, filter_val, val_cols, val_min,
        filter_params, params_cols, params_conditions)

    # Add Hue #
    baseline_params = {'data':df_data,  'x':x_column}
    baseline_params.update(resolve_hue(df_data, include_hue, hue_col, pal))

    # Autoannotate the Hits #
    if annot:
        assert isinstance(annot_label, str), "check param: autoannot_label"
        assert annot_label in df_data.columns.tolist(), "check param: autoannot_label"
        # annot_abs is a documented option but was missing from this gate, so
        # annot=True with only annot_abs set always tripped the assertion
        autoannot_bool = (isinstance(annot_top, int)
                          or isinstance(annot_cutoff, (int, float))
                          or isinstance(annot_abs, int))
        assert autoannot_bool, "check param: annot_top, annot_cutoff or annot_abs"

    fig, axes = plt.subplots(
        nrows=len(comparisons), ncols=1,
        **subplots_kws) # SETUP SUBPLOTS #
    if len(comparisons) == 1: axes = [axes]

    for ax, comp in zip(axes, comparisons):
        # Plot for Every Comparison #
        sns.scatterplot(ax=ax, y=comp, **baseline_params, **scatterplot_kws)

        # Overlay a Highlighted Subset #
        hits, _ = split_highlight(df_data, highlight, highlight_col, highlight_vals)
        if hits is not None and not hits.empty:
            ax.scatter(hits[x_column], hits[comp], **highlight_kws)

        # Overlay +/- n sd as Black Dashed Lines #
        draw_control_lines(ax, control_stdev(list_ctrlstats, comp),
                           n_sd=n_sd, **axhline_kws)

        # Annotate Hits #
        if annot:
            # select_annotated applies the same cutoff / top-n / abs-top-n
            # rules used by every other be_scan.plot function
            rows = select_annotated(df_data, comp, annot_top, annot_cutoff, annot_abs)
            if rows.empty and annot_top is None and annot_cutoff is None and annot_abs is None:
                warnings.warn("Please include a cutoff or top value")
            offset = 0.0 if annot_cutoff is not None else annot_offset
            annotate_points(ax, rows, x_column, comp, annot_label,
                            offset=offset, **label_text_kws)
            if annot_cutoff is None and not rows.empty:
                _display(rows)

        # Set Title and Axes #
        ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        style_axes(ax, xlab=xlab, ylab=ylab,
                   title=comp if title is None else title,
                   xlim_kws=xlim_kws, ylim_kws=ylim_kws,
                   hide_spines=('top', 'bottom', 'right'))
        draw_domains(ax, domains, domains_alpha, domains_color)

    # Save File #
    save_and_show(fig, savefig=savefig, show=show, out_dir=out_dir,
                  out_name=out_name, out_type=out_type, dpi=dpi)
    return fig, axes

def interactive_scatter(
    df_filepath, # DF Filename or DF #
    comparisons, # A List of Y-Axis Values #
    x_column, # X-Axis #

    annot_label='', # Label to Annotate Points #
    figsize=(600, 300),
    title="Scatter Plots", xlab='Amino Acid Position', ylab='sgRNA Score', # Scatterplot Labels #
    xwindow=None, # Filter Data on xwindow #
    include_hue=False, hue_col='', pal='pastel', # Color Parameters #
    domains=[], domains_alpha=0.25, domains_color='lightblue', # Draw Domains #
    neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # Negative Ctrl Parameters #
    savefig=True, show=True, out_name='scatterplot', out_dir='', # Output Parameters #

    # Style Params #
    xlim=None, ylim=None,
    plotly_go_kws={'opacity':0.85, 'line':dict(width=0.8)},

    # Added in 3.0.0 for parity with the other interactive plots. #
    filter_val=False, val_cols=[], val_min=0.0, # Filter Out Values #
    filter_params=False, params_cols=[], params_conditions=[['Missense', 'Silent', 'Mixed', 'Nonsense']], # Filter Out Categories #
    highlight=False, highlight_col='', highlight_vals=[], # Highlight a Subset #
    ):

    """[Summary]
    This function takes in a dataframe from analysis and then
    plots the data for each condition to reveal which guides are enriched

    Parameters
    ------------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    x_column : str, required
        column of .csv, typically amino acid position
    comparisons : list of str, required
        list of comparisons that correspond to columns of data

    annot_label : str, optional, defaults to ''
        name of an input column on which to annotate the interactive scatterplots
    figsize : tuple, optional, defaults to (1200, 600)
        (width, height) of figure in pixels
    title : str, optional, defaults to 'Scatter Plots'
        title of figure
    xlab : str, optional, defaults to 'Amino Acid Position'
        name of x-axis label
    ylab : str, optional, defaults to 'sgRNA Score'
        name of y-axis label
    xwindow : list, optional, defaults to None
        an inclusive range to restrict the x axis value to display

    include_hue: bool, optional, default to False
        whether or not to color points by a variable,
        will also restrict points plotted to only the hue_order values listed
    hue_col: str, optional, defaults to 'CtoT_muttype'
        the categorial dimension of the data, name of .csv data column
    pal: list of str, optional, defaults to a preset list of colors from ColorBrewer2
        a list of colors which correspond to hue_order

    domains: list, optional, defaults to []
        a list of domains with each item as {'start':n, 'end':m}
    domains_alpha: float, optional, defaults to 0.25
        alpha value of domains
    domains_color: str, optional, defaults to 'lightblue'
        color of domains

    neg_ctrl : bool, optional, defaults to False
        whether or not to calulate negative control for normalization and line drawing
    neg_ctrl_col : str, optional, defaults to 'CtoT_win_overlap'
        column of .csv which correspond to normalization control
    neg_ctrl_conditions : list of str, optional, defaults to ['Intron', 'Control']
        name of categories of neg_ctrl_col to normalize dataframe

    savefig : bool, optional, defaults to True
        whether or not to save the figure
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
    df_data = apply_xwindow(df_data, x_column, xwindow)
    df_data = apply_filters(
        df_data, filter_val, val_cols, val_min,
        filter_params, params_cols, params_conditions)

    # Add Hue #
    palette_dict = hue_color_map(df_data, include_hue, hue_col, pal)

    fig = make_subplots(
        rows=len(comparisons), cols=1,
        subplot_titles=comparisons, shared_xaxes=True)

    # Add Labels #
    if not annot_label in df_data.columns:
        warnings.warn(f"{annot_label} is not a column in input DataFrame", UserWarning)
        hovertext = None
    else: hovertext = df_data[annot_label]

    for idx, comp in enumerate(comparisons):
        # Scatter Plot #
        fig.add_trace(
            go.Scatter(
                x=df_data[x_column], y=df_data[comp],
                hovertext=hovertext, mode='markers', name=comp,
                marker=dict(size=6, color=df_data[hue_col].map(palette_dict) if include_hue else None),
                **plotly_go_kws),
            row=idx+1, col=1 )

        # Overlay +/- 2 sd as Black Dashed Lines #
        tup_comp_stdev = control_stdev(list_ctrlstats, comp)
        if tup_comp_stdev is not None:
            fig.add_shape(
                dict(type="line", x0=df_data[x_column].min(), x1=df_data[x_column].max(),
                     y0=2*tup_comp_stdev, y1=2*tup_comp_stdev,
                     line=dict(color="black", dash="dash")),
                row=idx+1, col=1 )
            fig.add_shape(
                dict(type="line", x0=df_data[x_column].min(), x1=df_data[x_column].max(),
                     y0=-2*tup_comp_stdev, y1=-2*tup_comp_stdev,
                     line=dict(color="black", dash="dash")),
                row=idx+1, col=1 )

        # Adjust X and Y Limits #
        x_range = [df_data[x_column].min()-10, df_data[x_column].max()+10] if xlim is None else [xlim['left'], xlim['right']]
        fig.update_xaxes(range=x_range, row=idx + 1, col=1 )
        if ylim is None:
            bound = max(abs(np.floor(df_data[comp].min())), abs(np.ceil(df_data[comp].max())))
            fig.update_yaxes(range=[-bound, bound], row=idx + 1, col=1 )
        else: fig.update_yaxes(range=ylim, row=idx + 1, col=1 )

        # Add Domain Highlighting #
        for d in domains:
            fig.add_shape(
                dict(type="rect", x0=d['start'], x1=d['end'],
                     y0=ylim[0] if ylim else -bound, y1=ylim[1] if ylim else bound,
                     fillcolor=domains_color, opacity=domains_alpha, layer="below",
                     ),
                row=idx+1, col=1 )

    # Update Layout #
    fig.update_layout(
        title_text=title,
        xaxis_title=xlab, yaxis_title=ylab,
        width=figsize[0], height=figsize[1] * len(comparisons),
        )

    # Save File #
    save_and_show_interactive(fig, savefig=savefig, show=show,
                              out_dir=out_dir, out_name=out_name)
    return fig
