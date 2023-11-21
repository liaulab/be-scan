"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from ._annotating_ import list_muttypes, color_list

def plot_corr_scatterplot(df_filepath, 
                          condition1, condition2, 
                          hue_column, 
                          hue_order=list_muttypes, palette=color_list, 
                          x_name='cond1 score', y_name='cond2 score', 
                          xmin=None, xmax=None, ymin=None, ymax=None, 
                          out_directory='', out_name='correlation_scatterplot', out_type='pdf', 
                          figsize=(4.5, 4), 
                          alpha=0.8, linewidth=1, edgecolor='black', s=25,
                          ):
    '''
    Make scatter plot for correlation assessment.
    '''

    df_data = pd.read_csv(df_filepath)
    df_filtered = df_data.loc[df_data[hue_column].isin(hue_order)]
    
    # Make plot
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(data=df_filtered, 
                    ax=ax, 
                    x=condition1, y=condition2, 
                    hue=hue_column, hue_order=hue_order, palette=palette, 
                    alpha=alpha, linewidth=linewidth, edgecolor=edgecolor, s=s
                    )
    
    # Adjust x and y axis limits
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    # Set labels
    plt.xlabel(x_name) # set x-axis label
    plt.ylabel(y_name) # set y-axis label
    # Adjust dimensions
    plt.tight_layout()
    
    # Save to pdf
    out = out_directory + condition1 + condition2 + '_' + out_name + '.' + out_type
    print(out)
    plt.savefig(out, format='pdf')
    plt.close()

# python3 -m be_scan plot_corr_scatterplot -df '../../../Downloads/NZL10196_v9_comparisons.csv' -c1 'd3-neg' -c2 'd9-pos' -hue 'Mut_type'
