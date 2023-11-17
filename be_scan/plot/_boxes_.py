"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import numpy as np
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt


def plot_boxes(df_input, cat_col, plot_x_list, y_val, # x y values, df, comparisons
               hue_order, palette, # sns.boxplot params
               plot_name, out_prefix, output_type='pdf', # output pdf params
               ylab='Log2 Fold Change', xlab='', # x and y label
               plot_col='comparison', # df params
               cutoff=20, swarm_list=None, jitter=True, # small categories and stripplot params
               dimensions=(5,4), list_negctrlstats=None, yax_set=True, # adjusting plots params
               swarm=False, # swarmplot
               sat_level=1, fliersize=4, width=0.4, # sns.boxplot params
               flierprops={'marker':'o', 'mec':'black', 'lw':1, 'alpha':0.8} # sns.boxplot params
               ):
    '''
    Make box plots.
    '''
    
    figpdf = PdfPages('_'.join([out_prefix, plot_name, 'boxes.pdf']))

    for comp in plot_x_list:
        
        # Isolate data corresponding to appropriate comparison/sample
        df_plot = df_input.loc[df_input[plot_col] == comp].copy()
        
        # Make boxplot
        fig, ax = plt.subplots(figsize=dimensions)
        sns.boxplot(x=cat_col, y=y_val, data=df_plot, 
                    order=hue_order, width=width, palette=palette, ax=ax, 
                    saturation=sat_level, fliersize=fliersize, flierprops=flierprops)
        plt.setp(ax.artists, edgecolor='black')
        plt.setp(ax.lines, color='black')

        # Check for small categories
        small_cats = []
        for category in hue_order:
            if df_plot[cat_col].value_counts()[category] < cutoff:
                small_cats.append(category)
        # If # of observations is below cutoff, stripplot instead of boxplot for those categories
        if len(small_cats) > 0:
            df_plot_s = df_plot.loc[df_plot[cat_col].isin(small_cats)].copy()
            df_plot = df_plot.loc[~df_plot[cat_col].isin(small_cats)]
            sns.stripplot(x=cat_col, y=y_val, data=df_plot_s, alpha=0.8,
                          order=hue_order, palette=palette, edgecolor='black',
                          jitter=jitter, linewidth=1, size=4, ax=ax)
            
        # make swarmplot if that is what user wants
        if swarm:
            if swarm_list != None:
                df_plot = df_plot.loc[df_plot[cat_col].isin(swarm_list)]
            sns.swarmplot(x=cat_col, y=y_val, data=df_plot,
                          order=hue_order, palette=palette, edgecolor='black',
                          linewidth=1, alpha=0.5, size=4, ax=ax)
            
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        if list_negctrlstats != None:
            tup_plot = [tup for tup in list_negctrlstats if tup[0] == comp][0]
            plt.axhline(y=2*tup_plot[2], color='k', ls='--', lw=1)
            plt.axhline(y=-2*tup_plot[2], color='k', ls='--', lw=1)
        
        # Adjust x and y axis limits
        if yax_set: 
            plt.ylim(np.floor(df_input[y_val].min()),
                     np.ceil(df_input[y_val].max()))
        
        # Set/adjust labels
        plt.title(comp) # Set plot title
        plt.ylabel(ylab) # Set y-axis label
        plt.xlabel(xlab) # Remove x-axis label
        plt.xticks(rotation=45, horizontalalignment='right')
        
        # Add major gridlines
        #plt.grid(b=True, which='major', axis='y', lw=0.5, c='k', alpha=0.2)
        # Adjust dimensions
        plt.tight_layout()
        
        # Save to pdf
        plt.savefig(figpdf, format=output_type)
        plt.close()

    figpdf.close()    
