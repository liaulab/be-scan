"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_corr_heatmap(df_input, plot_type, line_pos, plot_name, out_prefix,
                      cmap='coolwarm', center=None,
                      sns_context='paper', sns_palette='deep'):
    """
    Plots correlation heatmap and saves to pdf file.
    """
    # Plotting parameters and variables
    sns.set_context(sns_context)
    sns.set_palette(sns_palette)
    sns.set_style('ticks')
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Compute correlation matrix
    df_corr = df_input.corr(method=plot_type)

    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(4,4))
    ax = sns.heatmap(df_corr, cmap=cmap, center=center,
                     square=True, linewidths=0.5, cbar_kws={"shrink": 0.5}) #Heatmap

    # Rotate axis labels
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.yticks(rotation=0, horizontalalignment='right')
    for pos in line_pos:
        plt.hlines(pos, *ax.get_xlim())
        plt.vlines(pos, *ax.get_ylim())
    
    # Show frame
    for _, spine in ax.spines.items():
        spine.set_visible(True)
    
    # Adjustments
    plt.ylabel('')
    plt.xlabel('')
    plt.title(plot_name+' Correlation Heatmap')
    plt.tight_layout()

    # Save pdf and close everything
    plt.savefig('_'.join([out_prefix, plot_name, 'heatmap.pdf']), format='pdf')
    plt.close()
