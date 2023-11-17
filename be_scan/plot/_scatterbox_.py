"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot_scatterbox(df_scatter, hue_col, hue_order_s, palette_s,
                    df_boxes, cutoff, cat_col, hue_order_b, palette_b, jitter,
                    plot_list, list_negctrlstats, plot_name, out_prefix,
                    swarm=False):
    '''
    Make a scatterplot next to a boxplot with the y-axes aligned.

    '''
    
    #Plot
    figpdf = PdfPages('_'.join([out_prefix, plot_name, 'scatterbox.pdf']))

    for comp in plot_list:
        
        # Isolate data corresponding to appropriate comparison/sample
        df_plot1 = df_scatter.loc[df_scatter['comparison'] == comp].copy()
        df_plot2 = df_boxes.loc[df_boxes['comparison'] == comp].copy()
        
        # If number of observations is below cutoff, stripplot instead of
        # boxplot for those categories.
        small_cats = []
        for category in hue_order_b:
            if df_plot2[cat_col].value_counts()[category] < cutoff:
                small_cats.append(category)
        if len(small_cats) > 0:
            df_plot_s = df_plot2.loc[df_plot2[cat_col].isin(small_cats)].copy()
            df_plot2 = df_plot2.loc[~df_plot2[cat_col].isin(small_cats)]
        
        # Make plots
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(8,4),
                                       gridspec_kw={'width_ratios': [1.5, 1]})
        # Scatterplot on ax1
        sns.scatterplot(x='Edit_site_3A1', y='log2_fc', data=df_plot1, ax=ax1,
                        hue=hue_col, hue_order=hue_order_s, palette=palette_s,
                        alpha=0.8, linewidth=1, edgecolor='black', s=25)
        handles, labels = ax1.get_legend_handles_labels()
        legend = ax1.legend(handles=handles, labels=labels,
                            loc='lower left', bbox_to_anchor=(0.0, -0.02, 1, 0.3),
                            frameon=False, ncol=4)
        # Boxplot on ax2
        sns.boxplot(x=cat_col, y='log2_fc', data=df_plot2, order=hue_order_b,
                    palette=palette_b, saturation=1, fliersize=4, ax=ax2,
                    flierprops={'marker':'o', 'mec':'black', 'lw':1, 'alpha':0.8})
        plt.setp(ax2.artists, edgecolor='black')
        plt.setp(ax2.lines, color='black')
        if len(small_cats) > 0:
            sns.stripplot(x=cat_col, y='log2_fc', data=df_plot_s, alpha=0.8,
                          order=hue_order_b, palette=palette_b, jitter=jitter,
                          edgecolor='black', linewidth=1, size=4, ax=ax2)
        if swarm:
            sns.stripplot(x=cat_col, y='log2_fc', data=df_plot2,
                          order=hue_order_b, color='black',
                          jitter=0.1, alpha=0.5, size=2, ax=ax2)
        
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        tup_plot = [tup for tup in list_negctrlstats if tup[0] == comp][0]
        ax1.axhline(y=2*tup_plot[2], color='k', ls='--', lw=1)
        ax1.axhline(y=-2*tup_plot[2], color='k', ls='--', lw=1)
        ax2.axhline(y=2*tup_plot[2], color='k', ls='--', lw=1)
        ax2.axhline(y=-2*tup_plot[2], color='k', ls='--', lw=1)
        
        # Adjust x and y axis limits
        ax1.set_xlim(200, df_plot1['Edit_site_3A1'].max()+10)
        ax1.set_ylim(np.floor(df_scatter['log2_fc'].min()),
                     np.ceil(df_scatter['log2_fc'].max()))
                
        # Set labels
        ax1.set_title(comp)
        ax2.set_title(comp)
        ax1.set_xlabel('Amino Acid Position') #Set scatterplot x-axis label
        ax2.set_xlabel('') #Remove boxplot x-axis label
        for tick in ax2.get_xticklabels(): # Rotate boxplot x-axis labels
            tick.set_rotation(45)
            tick.set_horizontalalignment('right')
        ax1.set_ylabel('sgRNA Score') # Set y-axis label
        ax2.set_ylabel('')
        
        # Adjust dimensions
        plt.tight_layout()
        
        # Save to pdf
        plt.savefig(figpdf, format='pdf')
        plt.close()
    figpdf.close()    
