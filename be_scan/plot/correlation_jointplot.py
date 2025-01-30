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
import matplotlib.colors as mcolors

def corr_jointplot(df_filepath, condition1, condition2, 

    include_hue=False, hue_col='CtoT_muttype', pal='pastel', # color params
    savefig=True, show=True, out_dir='', out_name='correlation_jointplot', out_type='png', 

    # style params
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
    out_directory : str, optional, defaults to ''
        directory path of the output plot

    jointplot_kws : dict, optional, defaults to 
        {'alpha':0.8, 'linewidth':1, 
        'edgecolor':'black', 's':25}
        input params for sns.jointplot() 
        https://seaborn.pydata.org/generated/seaborn.jointplot.html

    Returns
    ------------
    """
    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)

    # set plot params, if there is hue add params for plotting
    baseline_params = {'data':df_data, 'x':condition1, 'y':condition2}
    if include_hue: 
        unique_types_sorted = sorted(df_data[hue_col].unique())
        baseline_params.update({'hue':hue_col, 'hue_order':unique_types_sorted})
        sns.set_palette(pal, len(unique_types_sorted))
    
    mpl.rcParams.update({'font.size': 10}) # STYLE #
    sns.jointplot(**baseline_params, **jointplot_kws, ) # PLOT #
    plt.tight_layout()
    ax = plt.gca()
    ax.legend(loc='center left', bbox_to_anchor=(1.25, 0.5)) # LEGEND #

    # calculate R2 stats and output them
    m, b, r_value, p_value, std_err = stats.linregress(df_data[condition1], df_data[condition2])
    print("R: {0} (p-value {1})".format(r_value, p_value))
    print("R2: {0}".format(r_value**2))

    # save to pdf and close
    outpath = Path(out_dir)
    if savefig: 
        out_name = f'{condition1}{condition2}_{out_name}.{out_type}'
        plt.savefig(outpath / out_name, format=out_type, dpi=300)
    if show: plt.show()
    plt.close()

def interactive_corr_jointplot(df_filepath, condition1, condition2, 

    include_hue=False, hue_col='CtoT_muttype', pal=px.colors.qualitative.Pastel, # color params
    savefig=True, show=True, out_dir='', out_name='correlation_jointplot', out_type='png', 
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
    out_directory : str, optional, defaults to ''
        directory path of the output plot

    Returns
    ------------
    """
    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)

    hex_pal = [mcolors.to_hex(c) for c in pal] 
    # Scatter plot with histograms
    fig = px.scatter(
        df_data, x=condition1, y=condition2, 
        color=hue_col if include_hue else None,  # Hue equivalent
        color_discrete_sequence=hex_pal if include_hue else None,  # Custom palette
        marginal_x="histogram",  marginal_y="histogram",  # Marginal histogram
        trendline="ols"  # Add regression line
    )

    # Regression statistics
    m, b, r_value, p_value, std_err = stats.linregress(df_data[condition1], df_data[condition2])
    r2 = r_value**2
    fig.add_annotation(
        x=min(df_data[condition1]), y=max(df_data[condition2]), 
        text=f"R² = {r2:.3f}<br>P-value = {p_value:.3g}", 
        showarrow=False, font=dict(size=12, color="black"), align="left", xanchor="left", yanchor="top"
    )

    # Update layout
    fig.update_layout(
        title=f"{condition1} vs {condition2}",
        xaxis_title=condition1, yaxis_title=condition2,
        legend_title=hue_col if include_hue else None,
        font=dict(size=10), 
        width=600, height=600, 
    )

    # Save file
    if show: fig.show()
    if savefig:
        outpath = Path(out_dir)
        out_name = f"{condition1}{condition2}_{out_name}.{out_type}"
        fig.write_html(str(outpath / out_name))
        