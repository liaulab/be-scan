"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function plots a scatterplot showing correlation between two given conditions}
"""

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import scipy.stats as stats
import plotly.express as px

def corr_jointplot(
    df_filepath, # DF Filename or DF #
    condition1, condition2, # Conditions to be Plotted #

    filter_val=False, val_cols=[], val_min=0.0, # Filter Out Values #
    filter_params=False, params_cols=[], params_conditions=[['Missense', 'Silent', 'Mixed', 'Nonsense']], # Filter Out Categories #
    include_hue=False, hue_col='CtoT_muttype', pal='pastel', # Color Parameters #
    savefig=True, show=True, out_dir='', out_name='correlation_jointplot', out_type='png', # Output Parameters #

    # Style Params #
    jointplot_kws={'alpha':0.8, 'linewidth':1, 'edgecolor':'black', 's':25},
    ):

    """[Summary]
    This function takes in a dataframe from count_reads, and plots
    a heatmap showing correlation between all given comparison conditions

    Parameters
    ------------
    df_filepath: str, required
        filepath to .csv data generated from count_reads
    condition1: str, required
        comparison condition 1, name of .csv data column
    condition2: str, required
        comparison condition 2, name of .csv data column

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

    include_hue: bool, optional, default to False
        whether or not to color points by a variable,
        will also restrict points plotted to only the hue_order values listed
    hue_col: str, optional, defaults to 'CtoT_muttype'
        the categorial dimension of the data, name of .csv data column
    palette: list of str, optional, defaults to a preset list of colors from ColorBrewer2
        a list of colors which correspond to hue_order

    savefig : bool, optional, defaults to True
        whether or not to save the figure
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_name : str, optional, defaults to 'scatterplot'
        name of the output plot
    out_type : str, optional, defaults to 'pdf'
        type of the output plot
    out_dir : str, optional, defaults to ''
        directory path of the output plot

    jointplot_kws : dict, optional, defaults to
        {'alpha':0.8, 'linewidth':1,
        'edgecolor':'black', 's':25}
        input params for sns.jointplot()
        https://seaborn.pydata.org/generated/seaborn.jointplot.html

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

    # Set up Plot Params #
    baseline_params = {'data':df_data, 'x':condition1, 'y':condition2}
    if include_hue:
        unique_types_sorted = sorted(df_data[hue_col].unique())
        baseline_params.update({'hue':hue_col, 'hue_order':unique_types_sorted})
        sns.set_palette(pal, len(unique_types_sorted))

    # Plot #
    sns.jointplot(**baseline_params, **jointplot_kws)
    ax = plt.gca()
    ax.legend(loc='center left', bbox_to_anchor=(1.25, 0.5))

    # Calculate R2 Stats #
    m, b, r_value, p_value, std_err = stats.linregress(df_data[condition1], df_data[condition2])
    print("R: {0} (p-value {1})".format(r_value, p_value))
    print("R2: {0}".format(r_value**2))

    plt.tight_layout()
    # Save File #
    outpath = Path(out_dir)
    if savefig:
        out = f'{condition1}{condition2}_{out_name}.{out_type}'
        plt.savefig(outpath / out, format=out_type, dpi=300)
    if show: plt.show()
    plt.close()

def interactive_corr_jointplot(
    df_filepath, # DF Filename or DF #
    condition1, condition2, # Conditions to be Plotted #

    annot_label='',
    figsize=(800,800),
    include_hue=False, hue_col='CtoT_muttype', pal=px.colors.qualitative.Pastel, # Color Parameters #
    savefig=True, show=True, out_dir='', out_name='correlation_jointplot', # Output Parameters #
    ):

    """[Summary]
    This function takes in a dataframe from count_reads, and plots
    a heatmap showing correlation between all given comparison conditions

    Parameters
    ------------
    df_filepath: str, required
        filepath to .csv data generated from count_reads
    condition1: str, required
        comparison condition 1, name of .csv data column
    condition2: str, required
        comparison condition 2, name of .csv data column

    include_hue: bool, optional, default to False
        whether or not to color points by a variable,
        will also restrict points plotted to only the hue_order values listed
    hue_col: str, optional, defaults to 'CtoT_muttype'
        the categorial dimension of the data, name of .csv data column
    pal: list of str, optional, defaults to a preset list of colors from ColorBrewer2
        a list of colors which correspond to hue_order
    annot_label : str, optional, defaults to ''
        name of an input column on which to annotate the interactive scatterplots
    figsize : tuple, optional, defaults to (1200, 600)
        (width, height) of figure in pixels

    savefig : bool, optional, defaults to True
        whether or not to save the figure
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_name : str, optional, defaults to 'scatterplot'
        name of the output plot
    out_dir : str, optional, defaults to ''
        directory path of the output plot

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

    hex_pal = ["#{:02x}{:02x}{:02x}".format(*tuple(map(int, c[4:-1].split(",")))) for c in pal]
    # Scatter plot with histograms
    fig = px.scatter(
        df_data, x=condition1, y=condition2,
        color=hue_col if include_hue else None,
        color_discrete_sequence=hex_pal if include_hue else None,
        marginal_x="histogram",  marginal_y="histogram",
        trendline="ols",
        hover_data={annot_label: True} if annot_label in df_data.columns else None
    )

    # Calculate Stats #
    m, b, r_value, p_value, std_err = stats.linregress(df_data[condition1], df_data[condition2])
    r2 = r_value**2
    fig.add_annotation(
        x=min(df_data[condition1]), y=max(df_data[condition2]),
        text=f"R² = {r2:.3f}<br>P-value = {p_value:.3g}",
        showarrow=False, font=dict(size=12, color="black"), align="left", xanchor="left", yanchor="top",
    )

    # Update Layout #
    fig.update_layout(
        title=f"{condition1} vs {condition2}",
        xaxis_title=condition1, yaxis_title=condition2,
        legend_title=hue_col if include_hue else None,
        font=dict(size=10), width=figsize[0], height=figsize[1],
    )

    # Save File #
    outpath = Path(out_dir)
    if savefig:
        fig.write_html(outpath / f"{condition1}{condition2}_{out_name}.html")
    if show: fig.show()
