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

import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


def plot_scatterplot(df_filepath, # dataframe
                     x_column, y_column, 
                     plot_column, hue_column, # df params
                     comparisons, 
                     neg_ctrl_category,
                     window=None,
                     plot_name='scatterplot', plot_type='pdf', out_directory='', # output params
                     xlab='Amino Acid Position', ylab='sgRNA Score', # scatterplot labels
                     alpha=0.8, linewidth=1, edgecolor='black', s=25, # scatterplot params
                     dimensions=(8,4)
                     ):
    """[Summary]
    :param [df_filepath]: [ParamDescription], defaults to [DefaultParamVal]
    :type [ParamName]: [ParamType](, optional)
    ...
    :raises [ErrorType]: [ErrorDescription]
    ...
    :return: [ReturnDescription]
    :rtype: [ReturnType]
    """

    in_dataframe = pd.read_csv(df_filepath)

    # Define list of annotations
    list_annocols = ['sgRNA_ID', 'Gene', 'Targeted_exon', 'C_count', 'is_C',
                    'Splice_check', 'Mut_type', 'Edit_site_3A1', 'Domain']
    # Normalize data to intergenic controls
    # calculate negative control stats
    df_negctrl, list_negctrlstats, avg_dict = calc_negative_controls(in_dataframe, comparisons, neg_ctrl_category)
    # perform normalization
    for comp in comparisons:
        in_dataframe[comp] = in_dataframe[comp].sub(avg_dict[comp])
    # tidy data
    list_cols = list_annocols + comparisons
    df_logfc = in_dataframe[list_cols].copy()
    for comp in comparisons: 
        df_logfc[comp+'_'+y_column] = in_dataframe[comp]

    # ColorBrewer2, 9 data classes, qualitative, 4th color scheme hex codes
    color_list = ['#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', 
                '#d9d9d9', '#8dd3c7', '#ffffb3', '#bebada', '#bc80bd']
    # Define lists to use for setting plotting parameters
    list_muttypes = ['Nonsense', 'Missense', 'Silent', 'Non-exon', 'Splice', 
                    'No_C/Exon', 'No_C/Non-exon', 'Control']
    hue_order_s = list_muttypes[:3]
    palette_s = color_list[:3]

    df_tidy = df_logfc.loc[(df_logfc[x_column]>=window[0]) & (df_logfc[x_column]<=window[1])]
    df_filtered = df_tidy.loc[df_tidy[hue_column].isin(list_muttypes[:3])]

    # output pdf information
    output_path = out_directory + plot_name + '.' + plot_type
    figpdf = PdfPages(output_path)
    
    for comp in comparisons:
        # Make plots
        fig, ax = plt.subplots(figsize=dimensions)
        
        # Scatterplot on ax1
        sns.scatterplot(ax=ax,
                        data=df_filtered, 
                        x=x_column, y=comp+'_'+y_column, 
                        hue=hue_column, 
                        hue_order=hue_order_s, palette=palette_s,
                        alpha=alpha, linewidth=linewidth, edgecolor=edgecolor, s=s)
        
        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        tup_plot = [tup for tup in list_negctrlstats if tup[0] == comp][0]
        ax.axhline(y=2*tup_plot[2], color=edgecolor, ls='--', lw=linewidth)
        ax.axhline(y=-2*tup_plot[2], color=edgecolor, ls='--', lw=linewidth)
        
        # Adjust x and y axis limits
        ax.set_xlim(df_filtered[x_column].min()-10, df_filtered[x_column].max()+10)
        ax.set_ylim(np.floor(df_filtered[comp+'_'+y_column].min()*1.5),
                     np.ceil(df_filtered[comp+'_'+y_column].max()*1.5))
        # Set tital and axes labels
        ax.set_title(comp)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)

        # Adjust dimensions
        plt.tight_layout()
        # Save to pdf
        plt.savefig(figpdf, format='pdf')
        plt.close()

    figpdf.close()  
    
### NORMALIZATION

# calculate the negative controls (ie the mean and stdev for the non)
def calc_negative_controls(df_data, list_compnames, neg_ctrl_category): 

    # Use negative controls to set cutoffs
    df_negctrl = df_data.loc[df_data['Gene'] == neg_ctrl_category].copy()
    list_negctrlstats = [] # list of tups of (comp, avg, sd, avg+2sd, avg-2sd)
    # cutoff_dict = {} # dictionary of comp: (avg+2sd, avg-2sd)
    avg_dict = {} # dictionary of comp: avg

    for comp in list_compnames:
        # mean and stdev
        temp = (df_negctrl[comp].mean(), df_negctrl[comp].std())
        # comparison, mean, 2 stdev above, 2 stdev below
        tup_comp = (comp, temp[0], temp[1], temp[0] + (2*temp[1]), temp[0] - (2*temp[1]))
        list_negctrlstats.append(tup_comp)
        # cutoff_dict[comp] = (tup_comp[3], tup_comp[4])
        avg_dict[comp] = temp[0]
        
    return df_negctrl, list_negctrlstats, avg_dict
