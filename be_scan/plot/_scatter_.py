"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196-3D-Clustering-v7.py Created on Mon Jun  1 19:03:10 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def plot_scatter(df_input, plot_list, hue_col, hue_order, palette,
                 list_negctrlstats, plot_name, out_prefix,
                 hue_norm=None, sns_context='paper', sns_palette='deep'):
    '''
    Make scatter plot.

    '''
    # Plotting parameters and variables
    sns.set_context(sns_context)
    sns.set_palette(sns_palette)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    #Plot
    figpdf = PdfPages('_'.join([out_prefix, plot_name, 'scatter.pdf']))

    for comp in plot_list:
        
        # Isolate data corresponding to appropriate comparison/sample
        df_plot = df_input.loc[df_input['comparison'] == comp].copy()
        
        # Make plot
        fig, ax = plt.subplots(figsize=(6,4))
        sns.scatterplot(x='Edit_site_3A1', y='log2_fc', data=df_plot, ax=ax,
                        hue=hue_col, hue_order=hue_order, palette=palette,
                        hue_norm=hue_norm, alpha=0.7, linewidth=1,
                        edgecolor='black', s=50)
        handles, labels = ax.get_legend_handles_labels()
        legend = plt.legend(handles=handles[1:], labels=labels[1:],
                            loc='lower right', bbox_to_anchor=(1.0, -0.02),
                            frameon=False, ncol=3)
        
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        tup_plot = [tup for tup in list_negctrlstats if tup[0] == comp][0]
        plt.axhline(y=tup_plot[3], color='k', ls='--', lw=1)
        plt.axhline(y=tup_plot[4], color='k', ls='--', lw=1)
        
        # Adjust x and y axis limits
        plt.xlim(df_plot['Edit_site_3A1'].min()-10, df_plot['Edit_site_3A1'].max()+10)
        plt.ylim(np.floor(df_input['log2_fc'].min()),
                 np.ceil(df_input['log2_fc'].max()))
        
        # Set labels
        plt.title(comp) # set plot title
        plt.xlabel('Amino Acid Position') # set x-axis label
        plt.ylabel('Log2 Fold Change') # set y-axis label
        # Add major gridlines
        plt.grid(b=True, which='major', axis='y', lw=0.5, c='k', alpha=0.2)
        # Adjust dimensions
        plt.tight_layout()
        
        # Save to pdf
        plt.savefig(figpdf, format='pdf')
        plt.close()
    figpdf.close()
