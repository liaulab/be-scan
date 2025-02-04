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

def scatterplot(df_filepath, # dataframe
                comparisons, # each comparison is a plot, and also the y axis
                x_column, # the x axis values

    xlab='Amino Acid Position', ylab='sgRNA Score', # scatterplot labels
    include_hue=False, hue_col='', pal='pastel', # color params
    neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # neg control params
    savefig=True, show=True, out_name='scatterplot', out_type='png', out_dir='', # output params
    domains=[], domains_alpha=0.25, domains_color='lightblue', # draw domains
    xwindow=[], 
    
    # style params
    subplots_kws={}, xlim={}, ylim={}, 
    xlim_kws={'xmin':None, 'xmax':None}, ylim_kws={'ymin':None, 'ymax':None},
    scatterplot_kws={'alpha': 0.85, 'linewidth': 0.8, 's': 30},
    axhline_kws={'color':'k', 'ls':'--', 'lw':1},
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
    include_hue: bool, optional, default to False
        whether or not to color points by a variable, 
        will also restrict points plotted to only the hue_order values listed
    hue_col: str, optional, defaults to 'CtoT_muttype'
        the categorial dimension of the data, name of .csv data column
    palette: list of str, optional, defaults to a preset list of colors from ColorBrewer2
        a list of colors which correspond to hue_order
        
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
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
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
        {'figsize':(4,4)}
        input params for plt.subplots() 
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
    axhline_kws : dict, optional, defaults to 
        {'color':'k', 'ls':'--', 'lw':1}
        input params for plt.axhline() 
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.axhline.html

    Returns
    ------------
    """
    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)

    # normalize data to intergenic controls if neg_ctrl is provided
    if neg_ctrl: 
        assert isinstance(neg_ctrl_col, str) and neg_ctrl_col in df_data.columns.tolist(), "check param: params_cols"
        assert isinstance(neg_ctrl_conditions, list), "check param: params_conditions"
        # calculate negative control stats
        _, list_ctrlstats, avg_dict = calc_neg_ctrls(df_data, comparisons, neg_ctrl_col, neg_ctrl_conditions)
        # calculate normalized log_fc scores for each comp condition
        df_data = norm_to_intergenic_ctrls(df_data, comparisons, avg_dict)

    if len(xwindow) == 2: 
        df_data = df_data[(df_data[x_column] >= xwindow[0]) & (df_data[x_column] <= xwindow[1])] # FILTER X WINDOW #
    elif len(xwindow) != 0: 
        warnings.warn(f"Warning: xwindow should be a list of 2 values", UserWarning)

    # if there is hue add params for plotting
    baseline_params = {'data':df_data,  'x':x_column}
    if include_hue: 
        unique_types_sorted = sorted(df_data[hue_col].unique())
        baseline_params.update({'hue':hue_col, 'hue_order':unique_types_sorted})
        sns.set_palette(pal, len(unique_types_sorted))

    mpl.rcParams.update({'font.size': 10}) # STYLE #
    fig, axes = plt.subplots(nrows=len(comparisons), ncols=1, figsize=(10, 3*len(comparisons)), 
                            **subplots_kws) # SETUP SUBPLOTS #
    if len(comparisons) == 1: axes = [axes]

    for ax, comp in zip(axes, comparisons):
        # make plots for every comparison
        sns.scatterplot(ax=ax, y=comp, **baseline_params, **scatterplot_kws)
        
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        if neg_ctrl and list_ctrlstats != None:
            tup_comp_stdev = [tup for tup in list_ctrlstats if tup[0] == comp][0][2]
            ax.axhline(y=2*tup_comp_stdev, **axhline_kws) # top baseline
            ax.axhline(y=-2*tup_comp_stdev, **axhline_kws) # bottom baseline

        # Adjust x and y axis limits
        if xlim == {}: ax.set_xlim(df_data[x_column].min()-10, df_data[x_column].max()+10)
        else: ax.set_xlim(**xlim)
        if ylim == {}: 
            bound = max(abs(np.floor(df_data[comp].min())), abs(np.ceil(df_data[comp].max())))
            ax.set_ylim(-1*bound, bound)
        else: ax.set_ylim(**ylim)
        # Set title and axes labels
        ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        ax.set_title(comp) ; ax.set_xlabel(xlab) ; ax.set_ylabel(ylab)
        ax.set_xlim(**xlim_kws) ; ax.set_ylim(**ylim_kws)
        ax.spines['top'].set_visible(False) ; ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        for d in domains:
            ax.axvspan(d['start'], d['end'], alpha=domains_alpha, facecolor=domains_color)

    plt.tight_layout() # ADJUST DIMENSIONS #
    # Save file
    outpath = Path(out_dir)
    if savefig: 
        out = f'{out_name}.{out_type}'
        plt.savefig(outpath / out, format=out_type, dpi=300)
    if show: plt.show()
    plt.close()

def interactive_scatter(df_filepath, # dataframe
                        comparisons, # each comparison is a plot, and also the y axis
                        x_column, # the x axis values

    xlab='Amino Acid Position', ylab='sgRNA Score', # scatterplot labels
    include_hue=False, hue_col='', pal='pastel', # color params
    neg_ctrl=False, neg_ctrl_col='', neg_ctrl_conditions=[], # neg control params
    savefig=True, show=True, out_name='scatterplot', out_dir='', # output params
    domains=[], domains_alpha=0.25, domains_color='lightblue', # draw domains
    xwindow=[], annot_label='', 
    
    # style params
    xlim={}, ylim={},
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

    xlab : str, optional, defaults to 'Amino Acid Position'
        name of x-axis label
    ylab : str, optional, defaults to 'sgRNA Score'
        name of y-axis label
    include_hue: bool, optional, default to False
        whether or not to color points by a variable, 
        will also restrict points plotted to only the hue_order values listed
    hue_col: str, optional, defaults to 'CtoT_muttype'
        the categorial dimension of the data, name of .csv data column
    palette: list of str, optional, defaults to a preset list of colors from ColorBrewer2
        a list of colors which correspond to hue_order
        
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
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
        path to output directory

    Returns
    ------------
    """
    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)

    # if there is hue add params for plotting
    if include_hue: 
        unique_types_sorted = sorted(df_data[hue_col].unique())
        sns.set_palette(pal, len(unique_types_sorted))

        palette = sns.color_palette(pal, n_colors=len(unique_types_sorted))
        pal = dict(zip(unique_types_sorted, [mcolors.to_hex(c) for c in palette]))

    # normalize data to intergenic controls if neg_ctrl is provided
    if neg_ctrl: 
        assert isinstance(neg_ctrl_col, str) and neg_ctrl_col in df_data.columns.tolist(), "check param: params_cols"
        assert isinstance(neg_ctrl_conditions, list), "check param: params_conditions"
        # calculate negative control stats
        _, list_ctrlstats, avg_dict = calc_neg_ctrls(df_data, comparisons, neg_ctrl_col, neg_ctrl_conditions)
        # calculate normalized log_fc scores for each comp condition
        df_data = norm_to_intergenic_ctrls(df_data, comparisons, avg_dict)

    if len(xwindow) == 2: 
        df_data = df_data[(df_data[x_column] >= xwindow[0]) & (df_data[x_column] <= xwindow[1])] # FILTER X WINDOW #
    elif len(xwindow) != 0: 
        warnings.warn(f"Warning: xwindow should be a list of 2 values", UserWarning)

    num_plots = len(comparisons)
    fig = make_subplots(rows=num_plots, cols=1, subplot_titles=comparisons, shared_xaxes=True)

    if not annot_label in df_data.columns: 
        warnings.warn(f"{annot_label} is not a column in input DataFrame", UserWarning)
        hovertext = None
    else: hovertext = df_data[annot_label]
    
    # Iterate through comparisons to create each subplot
    for idx, comp in enumerate(comparisons): 
        # Scatter Plot
        fig.add_trace(
            go.Scatter(
                x=df_data[x_column], y=df_data[comp], 
                hovertext=hovertext, mode='markers', name=comp, 
                marker=dict(size=6, color=df_data[hue_col].map(pal) if include_hue else None), 
                **plotly_go_kws), row=idx + 1, col=1 )

        # Overlay neg ctrl avg +/- 2 SD as black dashed lines
        if neg_ctrl and list_ctrlstats != None:
            tup_comp_stdev = next(tup[2] for tup in list_ctrlstats if tup[0] == comp)
            fig.add_shape(
                dict(type="line", x0=df_data[x_column].min(), x1=df_data[x_column].max(),
                     y0=2*tup_comp_stdev, y1=2*tup_comp_stdev, line=dict(color="black", dash="dash")),
                row=idx + 1, col=1 )
            fig.add_shape(
                dict(type="line", x0=df_data[x_column].min(), x1=df_data[x_column].max(),
                     y0=-2*tup_comp_stdev, y1=-2*tup_comp_stdev, line=dict(color="black", dash="dash")),
                row=idx + 1, col=1 )

        # Adjust x and y axis limits
        x_range = [df_data[x_column].min()-10, df_data[x_column].max()+10] if not xlim else [xlim['left'], xlim['right']]
        fig.update_xaxes(range=x_range, row=idx + 1, col=1 )
        if not ylim:
            bound = max(abs(np.floor(df_data[comp].min())), abs(np.ceil(df_data[comp].max())))
            fig.update_yaxes(range=[-bound, bound], row=idx + 1, col=1 )
        else: 
            fig.update_yaxes(range=ylim, row=idx + 1, col=1 )
        
        # Add domain highlighting
        for d in domains:
            fig.add_shape(
                dict(type="rect", x0=d['start'], x1=d['end'], y0=ylim[0] if ylim else -bound, 
                     y1=ylim[1] if ylim else bound, fillcolor=domains_color, opacity=domains_alpha, layer="below"),
                row=idx + 1, col=1 )

    # Update layout
    fig.update_layout(
        title_text="Scatter Plots - Comparisons", 
        xaxis_title=xlab, yaxis_title=ylab, 
        width=1200, height=600 * num_plots, )

    # Save file
    if show: fig.show()
    if savefig:
        outpath = Path(out_dir)
        out_name = f"{out_name}"
        fig.write_html(str(outpath / out_name))

# scatterplot(
#     df_filepath="tests/test_data/plot/NZL10196_v9_comparisons.csv", 
#     comparisons=["d3-pos", "d3-neg"], 
#     x_column='Edit_site_3A1', 
#     include_hue=True, hue_col='Mut_type', 
#     neg_ctrl=True, neg_ctrl_col='Gene', neg_ctrl_conditions=['NON-GENE'], # neg control params
#     # xlim={'left':200, 'right':920}, 
#     # annot=True, annot_label='sgRNA_ID', annot_abs=10, # annot_cutoff=0.5 annot_top=10
#     # domains=[{'start':300, 'end':400}], 
#     xwindow=[400,700]
# )
# scatterplot(
#     df_filepath="tests/test_data/plot/NZL10196_v9_comparisons.csv", 
#     comparisons=["d3-pos", "d3-neg"], 
#     x_column='Edit_site_3A1', 
#     include_hue=True, hue_col='Mut_type', 
#     neg_ctrl=True, neg_ctrl_col='Gene', neg_ctrl_conditions=['NON-GENE'], # neg control params
#     xlim={'left':200, 'right':920}, 
#     # annot=True, annot_label='sgRNA_ID', annot_abs=10, # annot_cutoff=0.5 annot_top=10
#     # domains=[{'start':300, 'end':400}], 
#     # interactive=True, 
# )
