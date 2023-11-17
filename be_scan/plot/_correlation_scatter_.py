"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import seaborn as sns
import matplotlib.pyplot as plt

def plot_corr_scatter(df_input, x_col, x_name, y_col, y_name, out_prefix,
                      out_suffix, hue_col, hue_order, palette,
                      xmin=None, xmax=None, ymin=None, ymax=None
                      ):
    '''
    Make scatter plot for correlation assessment.
    '''
    
    # Make plot
    fig, ax = plt.subplots(figsize=(4.5,4))
    sns.scatterplot(x=x_col, y=y_col, data=df_input, ax=ax,
                    hue=hue_col, hue_order=hue_order, palette=palette,
                    alpha=0.7, linewidth=1, edgecolor='black', s=25)
    handles, labels = ax.get_legend_handles_labels()
    legend = plt.legend(handles=handles[1:], labels=labels[1:])
    
    # Adjust x and y axis limits
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    # Set labels
    plt.xlabel(x_name) # set x-axis label
    plt.ylabel(y_name) # set y-axis label
    # Add major gridlines
    #plt.grid(b=True, which='major', axis='y', lw=0.5, c='k', alpha=0.2)
    # Adjust dimensions
    plt.tight_layout()
    
    # Save to pdf
    plt.savefig('_'.join([out_prefix, x_col, y_col,
                          'correlation_by', out_suffix, '.pdf']), format='pdf')
    plt.close()
