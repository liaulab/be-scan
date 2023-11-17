"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import seaborn as sns
import matplotlib.pyplot as plt

def plot_corr_heatmap(df_input, plot_type, 
                      out_prefix, plot_name, xlab='', ylab='', # output related params
                      line_pos=[2,4], cmap='coolwarm', center=None, linewidth=0.5 # plot related params
                     ):
    """
    Plots correlation heatmap and saves to pdf file.
    """
    # Plotting parameters and variables
    sns.set_style('ticks')
    
    # Compute correlation matrix
    df_corr = df_input.corr(method=plot_type)

    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(4,4))
    ax = sns.heatmap(df_corr, cmap=cmap, center=center,
                     square=True, linewidths=linewidth, cbar_kws={"shrink": 0.5}) # heatmap

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
    plt.ylabel(xlab)
    plt.xlabel(ylab)
    plt.title(plot_name)
    plt.tight_layout()

    # Save pdf and close everything
    plt.savefig('_'.join([out_prefix, plot_name, 'heatmap.pdf']), format='pdf')
    plt.close()
