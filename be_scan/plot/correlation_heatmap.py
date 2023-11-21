"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

def plot_corr_heatmap(df_filepath, 
                      plot_type, 
                      comparisons, 
                      out_directory='', out_name='correlation_heatmap', out_type='pdf', 
                      xlab='', ylab='', title='Spearman Correlation Heatmap', # output related params
                      line_pos=[2,4], cmap='coolwarm', center=0, linewidth=0.5 # plot related params
                     ): 
    """
    Plots correlation heatmap and saves to pdf file.
    """

    # Plotting parameters and variables
    sns.set_style('ticks')

    # Compute correlation matrix
    df = pd.read_csv(df_filepath)
    df_comp = df[comparisons].copy()
    df_corr = df_comp.corr(method=plot_type)

    # Set up the matplotlib figure
    fig, ax = plt.subplots(figsize=(4,4))
    ax = sns.heatmap(df_corr, 
                     cmap=cmap, center=center,
                     square=True, linewidths=linewidth, cbar_kws={"shrink": 0.5}
                     )

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
    plt.title(title)

    plt.tight_layout()
    # Save pdf and close everything
    output_path = out_directory + out_name + '.' + out_type
    plt.savefig(output_path, format='pdf')
    plt.close()

# python3 -m be_scan plot_corr_heatmap -df '../../../Downloads/NZL10196_v9_comparisons.csv' -p 'spearman' -c 'd3-pos' 'd3-neg' 'd6-pos' 'd6-neg' 'd9-pos' 'd9-neg'
