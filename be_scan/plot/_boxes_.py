"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import numpy as np
import seaborn as sns
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def plot_boxes(df_input, plot_list, cutoff, cat_col, hue_order, palette,
               plot_name, out_prefix, width=0.8, size=(5,4),
               swarm=False, swarm_list=None,
               plot_col='comparison', list_negctrlstats=None, yax_set=True,
               y_col='log2_fc', y_label='Log2 Fold Change', jitter=True,
               sns_context='paper', sns_palette='deep'):
    '''
    Make box plots.

    '''
    # Plotting parameters and variables
    sns.set_context(sns_context)
    sns.set_palette(sns_palette)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    #Plot
    figpdf = PdfPages('_'.join([out_prefix, plot_name, 'boxes.pdf']))

    for comp in plot_list:
        
        # Isolate data corresponding to appropriate comparison/sample
        df_plot = df_input.loc[df_input[plot_col] == comp].copy()
        
        # If number of observations is below cutoff, stripplot instead of
        # boxplot for those categories.
        small_cats = []
        for category in hue_order:
            if df_plot[cat_col].value_counts()[category] < cutoff:
                small_cats.append(category)
        if len(small_cats) > 0:
            df_plot_s = df_plot.loc[df_plot[cat_col].isin(small_cats)].copy()
            df_plot = df_plot.loc[~df_plot[cat_col].isin(small_cats)]
        
        # Make plot
        fig, ax = plt.subplots(figsize=size)
        sns.boxplot(x=cat_col, y=y_col, data=df_plot, order=hue_order,
                    width=width, showfliers=False,
                    palette=palette, saturation=1, ax=ax)
        plt.setp(ax.artists, edgecolor='black')
        plt.setp(ax.lines, color='black')
        if len(small_cats) > 0:
            sns.stripplot(x=cat_col, y=y_col, data=df_plot_s, alpha=0.8,
                          order=hue_order, palette=palette, edgecolor='black',
                          jitter=jitter, linewidth=1, size=4, ax=ax)
        if swarm:
            if swarm_list != None:
                df_plot = df_plot.loc[df_plot[cat_col].isin(swarm_list)]
            sns.stripplot(x=cat_col, y=y_col, data=df_plot,
                          order=hue_order, color='black', jitter=jitter,
                          alpha=0.5, size=4, ax=ax)
            
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        if list_negctrlstats != None:
            tup_plot = [tup for tup in list_negctrlstats if tup[0] == comp][0]
            plt.axhline(y=2*tup_plot[2], color='k', ls='--', lw=1)
            plt.axhline(y=-2*tup_plot[2], color='k', ls='--', lw=1)
        
        # Adjust x and y axis limits
        if yax_set:
            plt.ylim(np.floor(df_input[y_col].min()),
                     np.ceil(df_input[y_col].max()))
        
        # Set/adjust labels
        plt.title(comp) # Set plot title
        plt.ylabel(y_label) # Set x-axis label
        plt.xlabel('') # Remove y-axis label
        plt.xticks(rotation=45, horizontalalignment='right')
        
        # Add major gridlines
        #plt.grid(b=True, which='major', axis='y', lw=0.5, c='k', alpha=0.2)
        # Adjust dimensions
        plt.tight_layout()
        
        # Save to pdf
        plt.savefig(figpdf, format='pdf')
        plt.close()
    figpdf.close()
