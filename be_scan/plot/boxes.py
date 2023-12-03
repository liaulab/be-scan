"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: This function plots chosen guides by plot_column categories to show the distribution of guides}
"""

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

from be_scan.plot._annotating_ import norm_to_intergenic_ctrls, calc_negative_controls

def plot_boxes(df_filepath, 
               plot_column, plot_conditions, 
               y_column, 
               comparisons, 
               neg_ctrl_col, neg_ctrl_category,

               filter_column='Mut_type', filter_category='Missense',
               xlab='', ylab='Log2 Fold Change', # x and y label
               out_name='boxes', out_type='pdf', out_directory='', # output params
               dimensions=(5,4), 
               yax_set=True, # adjusting plots params
               sat_level=1, fliersize=4, width=0.4, # sns.boxplot params
               flierprops={'marker':'o', 'mec':'black', 'lw':1, 'alpha':0.8}, # sns.boxplot params
               savefig=True,
               ):
    
    """[Summary]
    This function takes in a dataframe from count_reads, performs normalization, 
    and then plots chosen (ex control) guides by plot_column categories, 
    to show the distribution of categories of guides
    
    Parameters
    ----------
    df_filepath : str, required
        filepath to .csv data generated from count_reads
    plot_columnstr, required
        column of .csv, typically domain or mutation type
    plot_conditions : list of str, required
        category names of plot_column
    y_column : str, required
        column of .csv, typically the normalized log_fc score
    comparisons : list of str, required
        list of comparisons that correspond to columns of data
    neg_ctrl_col : str, required
        column of .csv which correspond to normalization control
    neg_ctrl_category : str, required
        categorical variable of neg_ctrl_col

    filter_column : str, optional, defaults to 'Mut_type'
        name of column to filter dataframe for plotting
    filter_category : str, optional, defaults to 'Missense'
        name of categories of filter_column to filter dataframe
    xlab : str, optional, defaults to ''
        name of x-axis label
    ylab : str, optional, defaults to 'Log2 Fold Change'
        name of y-axis label
    out_name : str, optional, defaults to 'scatterplot'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
        path to output directory
    savefig : boolean, optional, defaults to True
        option of saving figure to output or not

    Returns
    ----------
    None
    """

    df_input = pd.read_csv(df_filepath)

    # Normalize data to intergenic controls
    # calculate negative control stats
    _, list_negctrlstats, avg_dict = calc_negative_controls(df_input, comparisons, neg_ctrl_col, neg_ctrl_category)
    # calculate normalized log_fc scores for each comp condition
    df_logfc = norm_to_intergenic_ctrls(df_input, comparisons, avg_dict, y_column)
    
    df_filtered = df_logfc.loc[df_logfc[filter_column]==filter_category]
    df_filtered = df_filtered.loc[df_filtered[plot_column].isin(plot_conditions)].copy()
    
    for comp in comparisons:
                
        # Make boxplot
        fig, ax = plt.subplots(figsize=dimensions)
        sns.boxplot(data=df_filtered, 
                    ax=ax, 
                    x=plot_column, y=comp+'_'+y_column, 
                    # order=hue_order, palette=palette, 
                    width=width, saturation=sat_level, fliersize=fliersize, flierprops=flierprops)
        plt.setp(ax.artists, edgecolor='black')
        plt.setp(ax.lines, color='black')

        # Overlay neg ctrl avg +/- 2 sd as black dashed line
        if list_negctrlstats != None:
            tup_plot = [tup for tup in list_negctrlstats if tup[0] == comp][0]
            plt.axhline(y=2*tup_plot[2], color='k', ls='--', lw=1)
            plt.axhline(y=-2*tup_plot[2], color='k', ls='--', lw=1)
        
        # Adjust x and y axis limits
        if yax_set: 
            plt.ylim(np.floor(df_filtered[comp+'_'+y_column].min()),
                     np.ceil(df_filtered[comp+'_'+y_column].max()))
        
        # Set/adjust labels
        plt.title(comp) # Set plot title
        plt.ylabel(ylab) # Set y-axis label
        plt.xlabel(xlab) # Remove x-axis label
        plt.xticks(rotation=45, horizontalalignment='right')
        
        # Adjust dimensions
        plt.tight_layout()
        # Save to pdf
        if savefig: 
            output_path = out_directory + out_name + comp + '.' + out_type
            plt.savefig(output_path, format=out_type)
        plt.show()
        plt.close()

# python3 -m be_scan plot_boxes -df '../../../Downloads/NZL10196_v9_comparisons.csv' -p 'Domain' -pc 'PWWP' 'ADD' 'MTase' -y 'log2_fc' -c 'd3-pos' 'd3-neg' 'd6-pos' 'd6-neg' 'd9-pos' 'd9-neg' -ncol 'Gene' -ncat "NON-GENE"
