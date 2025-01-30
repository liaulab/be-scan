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
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

def sgrna_enrichment(df_filepath, comparisons, 

    highlight=False, highlight_col='sgRNA_ID', highlight_vals=['sgRNA_0'], 
    density=True, 
    
    savefig=True, show=True, out_name='scatterplot', out_type='png', out_directory='', # output params
    interactive=False, 

    rugplot_kws={'height':0.25, 'color':'black', 'linewidth':2, 'alpha':0.05}, 
    kdeplot_kws={'color':'black', 'linewidth':1, 'alpha':0.5}, 
    highlight_rugplot_kws={'height':0.25, 'color':'red', 'linewidth':5, 'alpha':0.6}
    ):
    
    """[Summary]
    This function takes in a dataframe from analysis and then 
    plots the data for each condition to highlight specifically enriched guides

    Parameters
    ------------
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
    
    if interactive: 
        # Set up subplots if multiple comparisons exist
        num_comparisons = len(comparisons)
        fig = make_subplots(rows=1, cols=num_comparisons, subplot_titles=[f"sgRNA Enrichment for {val}" for val in comparisons])

        # Loop through each comparison to add density and rug plots
        for i, val in enumerate(comparisons):
            # Add density plot (Kernel Density Estimation)
            if density:
                kde_fig = px.histogram(df_data, x=val, marginal="density", histnorm='probability density', opacity=0.6)
                for trace in kde_fig.data:  # Extract traces from px.histogram and add to our figure
                    fig.add_trace(trace, row=1, col=i+1)
            
            # Add rug plot (small markers along the x-axis)
            fig.add_trace(go.Scatter(
                x=df_data[val], y=[-0.01] * len(df_data),  # Small offset below the density plot
                mode="markers", marker=dict(symbol="line-ns-open", color="black", size=6),
                name=f"Rugplot {val}",
            ), row=1, col=i+1)

            # Highlight specific guides if needed
            if highlight:
                highlight_df = df_data[df_data[highlight_col].isin(highlight_vals)]
                fig.add_trace(go.Scatter(
                    x=highlight_df[val], y=[-0.02] * len(highlight_df),  # Slightly lower offset for distinction
                    mode="markers",
                    marker=dict(symbol="line-ns-open", color="red", size=6),
                    name=f"Highlighted {val}",
                ), row=1, col=i+1)

        # Update layout settings
        fig.update_layout(
            title="sgRNA Enrichment", showlegend=True,
            height=400, width=800 * num_comparisons,
            xaxis_title="Enrichment Values", yaxis_title="Density",
        )
 
        # Show plot
        if show: fig.show()
        # Save the figure
        outpath = Path.cwd() / out_directory
        if savefig:
            out_name_full = f"{out_name}.html"
            fig.write_html(outpath / out_name_full)  # Saves as an interactive HTML

    else: 
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
    