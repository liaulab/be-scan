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
    if isinstance(df_filepath, pd.DataFrame):
        df_data = df_filepath.copy()
    elif isinstance(df_filepath, (str, Path)):
        df_filepath = Path(df_filepath)
        df_data = pd.read_csv(df_filepath)
    else:
        raise ValueError("df_filepath must be a pandas DataFrame or a path string/Path.")

    # Normalize Data to Intergenic Ctrls #
    if neg_ctrl:
        assert isinstance(neg_ctrl_col, str) and neg_ctrl_col in df_data.columns.tolist(), "check param: params_cols"
        assert isinstance(neg_ctrl_conditions, list), "check param: params_conditions"
        # Calculate Negative Ctrl Stats #
        _, list_ctrlstats, avg_dict = calc_neg_ctrls(df_data, comparisons, neg_ctrl_col, neg_ctrl_conditions)
        # Normalize Each Comparison #
        df_data = norm_to_intergenic_ctrls(df_data, comparisons, avg_dict)

    # Filter Data Before Plotting #
    if xwindow is not None and len(xwindow) == 2:
        df_data = df_data[(df_data[x_column] >= xwindow[0]) & (df_data[x_column] <= xwindow[1])]
    elif xwindow is not None and len(xwindow) != 2:
        warnings.warn(f"Warning: xwindow should be a list of 2 values", UserWarning)

    # Add Hue #
    baseline_params = {'data':df_data,  'x':x_column}
    if include_hue:
        assert isinstance(hue_col, str) and hue_col in df_data.columns, "check param: hue_col"
        unique_types = df_data[hue_col].dropna().unique()
        baseline_params.update({'hue': hue_col, 'hue_order': unique_types, 'palette': sns.color_palette(pal, len(unique_types))})

    # Filter by Filter Params #
    if filter_val:
        assert isinstance(val_min, float), "check param: val_min"
        assert isinstance(val_cols, list) and len(val_cols) > 0, "check param: val_cols"
        for col in val_cols:
            print(f'Filter data by {col} to be greater than {str(val_min)}')
            df_data = df_data[df_data[col] > val_min]
    if filter_params:
        assert isinstance(params_cols, list), "check param: params_cols"
        assert isinstance(params_conditions, list), "check param: params_conditions"
        for col, conds in zip(params_cols, params_conditions):
            print(f'Filter data by {col} to contain:')
            print(conds)
            df_data = df_data.loc[df_data[col].isin(conds)]

    # Autoannotate the Hits #
    if annot:
        assert isinstance(annot_label, str), "check param: autoannot_label"
        assert annot_label in df_data.columns.tolist(), "check param: autoannot_label"
        autoannot_bool = isinstance(annot_top, int) or isinstance(annot_cutoff, float)
        assert autoannot_bool, "check param: autoannot_top autoannot_cutoff"

    fig, axes = plt.subplots(
        nrows=len(comparisons), ncols=1,
        **subplots_kws) # SETUP SUBPLOTS #
    if len(comparisons) == 1: axes = [axes]

    for ax, comp in zip(axes, comparisons):
        # Plot for Every Comparison #
        sns.scatterplot(ax=ax, y=comp, **baseline_params, **scatterplot_kws)

        # Overlay +/- 2 sd as Black Dashed Lines #
        if neg_ctrl and list_ctrlstats != None:
            tup_comp_stdev = [tup for tup in list_ctrlstats if tup[0] == comp][0][2]
            ax.axhline(y=2*tup_comp_stdev, **axhline_kws)
            ax.axhline(y=-2*tup_comp_stdev, **axhline_kws)

        # Annotate Hits #
        if annot:
            if annot_cutoff:
                # Annotate All Above the Cutoff #
                df_cutoff = df_data[df_data[comp].abs() > annot_cutoff]
                for _, row in df_cutoff.iterrows():
                    if row[annot_label] != "nan":
                        ax.text(row[x_column], row[comp], row[annot_label], **label_text_kws)
            elif annot_abs:
                # Annotate the Top n Scoring #
                sorted_df = df_data.assign(yabs=df_data[comp].abs()).sort_values(by='yabs', ascending=False)
                toprows_df = sorted_df.head(annot_abs)
                for _, row in toprows_df.iterrows():
                    t = ax.text(row[x_column], row[comp]-0.2, row[annot_label], **label_text_kws)
                display(toprows_df)
            elif annot_top:
                # Annotate the Top or Bottom n #
                sign = np.sign(annot_top)
                ascen = True if sign==-1 else False
                sorted_df = df_data.assign(yabs=df_data[comp]).sort_values(by='yabs', ascending=ascen)
                toprows_df = sorted_df.head(sign*annot_top)
                for _, row in toprows_df.iterrows():
                    t = ax.text(row[x_column], row[comp]-0.2, row[annot_label], **label_text_kws)
                display(toprows_df)
            else: warnings.warn("Please include a cutoff or top value")

        # Set Title and Axes #
        ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        ax.set_title(comp) ; ax.set_xlabel(xlab) ; ax.set_ylabel(ylab)
        ax.set_xlim(**xlim_kws) ; ax.set_ylim(**ylim_kws)
        ax.spines['top'].set_visible(False) ; ax.spines['bottom'].set_visible(False) ; ax.spines['right'].set_visible(False)
        for d in domains:
            ax.axvspan(d['start'], d['end'], alpha=domains_alpha, facecolor=domains_color)

    plt.tight_layout()
    # Save File #
    outpath = Path(out_dir)
    if savefig:
        out = f'{out_name}.{out_type}'
        plt.savefig(outpath / out, format=out_type, dpi=300)
    if show: plt.show()
    plt.close()

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
    if isinstance(df_filepath, pd.DataFrame):
        df_data = df_filepath.copy()
    elif isinstance(df_filepath, (str, Path)):
        df_filepath = Path(df_filepath)
        df_data = pd.read_csv(df_filepath)
    else:
        raise ValueError("df_filepath must be a pandas DataFrame or a path string/Path.")

    # Normalize Data to Intergenic Ctrls #
    if neg_ctrl:
        assert isinstance(neg_ctrl_col, str) and neg_ctrl_col in df_data.columns.tolist(), "check param: params_cols"
        assert isinstance(neg_ctrl_conditions, list), "check param: params_conditions"
        # Calculate Negative Ctrl Stats #
        _, list_ctrlstats, avg_dict = calc_neg_ctrls(df_data, comparisons, neg_ctrl_col, neg_ctrl_conditions)
        # Normalize Each Comparison #
        df_data = norm_to_intergenic_ctrls(df_data, comparisons, avg_dict)

    # Filter Data Before Plotting #
    if xwindow is not None and len(xwindow) == 2:
        df_data = df_data[(df_data[x_column] >= xwindow[0]) & (df_data[x_column] <= xwindow[1])]
    elif xwindow is not None and len(xwindow) != 2:
        warnings.warn(f"Warning: xwindow should be a list of 2 values", UserWarning)

    # Add Hue #
    if include_hue:
        assert isinstance(hue_col, str) and hue_col in df_data.columns.tolist(), "check param: hue_col"
        unique_types = df_data[hue_col].dropna().unique()
        pal_list = sns.color_palette(pal, n_colors=len(unique_types))
        palette_dict = dict(zip(unique_types, [mcolors.to_hex(c) for c in pal_list]))

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
        if neg_ctrl and list_ctrlstats != None:
            tup_comp_stdev = next(tup[2] for tup in list_ctrlstats if tup[0] == comp)
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
    outpath = Path(out_dir)
    if savefig:
        fig.write_html(str(outpath / f"{out_name}.html"))
    if show: fig.show()
