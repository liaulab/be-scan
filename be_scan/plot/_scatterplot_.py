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

def plot_scatterplot(df_scatter, # dataframe
                     plot_col, hue_col, # df params
                     xval, yval, xlab, ylab, # scatterplot basic params
                     plot_x_list, list_negctrlstats, # lists of inputs
                     plot_name, out_prefix, # output params
                     hue_order_s, palette_s, alpha=0.8, linewidth=1, edgecolor='black', s=25, # scatterplot params
                     dimensions=(8,4)
                     ):
    '''
    Make a scatterplot
    '''
    
    #Plot
    figpdf = PdfPages('_'.join([out_prefix, plot_name, 'scatterbox.pdf']))

    for comp in plot_x_list:
        
        # Isolate data corresponding to appropriate comparison/sample
        df_plot1 = df_scatter.loc[df_scatter[plot_col] == comp].copy()
        
        # Make plots
        fig, ax = plt.subplots(figsize=dimensions)
        
        # Scatterplot on ax1
        sns.scatterplot(x=xval, y=yval, data=df_plot1, ax=ax,
                        hue=hue_col, hue_order=hue_order_s, palette=palette_s,
                        alpha=alpha, linewidth=linewidth, edgecolor=edgecolor, s=s)
        
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        tup_plot = [tup for tup in list_negctrlstats if tup[0] == comp][0]
        ax.axhline(y=2*tup_plot[2], color='k', ls='--', lw=1)
        ax.axhline(y=-2*tup_plot[2], color='k', ls='--', lw=1)
        
        # Adjust x and y axis limits
        ax.set_xlim(200, df_plot1[xval].max()+10)
        ax.set_ylim(np.floor(df_scatter[yval].min()),
                     np.ceil(df_scatter[yval].max()))
                
        # Set labels
        ax.set_title(comp)
        ax.set_xlabel(xlab) #Set scatterplot x-axis label
        ax.set_ylabel(ylab) # Set y-axis label
        
        # Adjust dimensions
        plt.tight_layout()
        
        # Save to pdf
        plt.savefig(figpdf, format='pdf')
        plt.close()

    figpdf.close()    
