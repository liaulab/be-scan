"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: some helper functions and helper reference variables relevant for plotting}
"""

# ColorBrewer2, 9 data classes, qualitative, 4th color scheme hex codes
color_list = ['#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', 
            '#d9d9d9', '#8dd3c7', '#ffffb3', '#bebada', '#bc80bd']
# Define lists to use for setting plotting parameters
list_muttypes = ['Nonsense', 'Missense', 'Silent', 'Non-exon', 'Splice', 
                'No_C/Exon', 'No_C/Non-exon', 'Control']

def norm_to_intergenic_ctrls(in_dataframe, comparisons, avg_dict): 
    """[Summary]
    This function normalizes data in in_dataframe to a 
    set of controls calculated from calc_negative_control, 
    typically the intergenic controls (control guides that do not target the gene)

    Parameters
    ----------

    Returns
    ----------
    df_logfc : pandas dataframe
    """
    # perform normalization for each cond individually
    for comp in comparisons:
        in_dataframe[comp] = in_dataframe[comp].sub(avg_dict[comp])

    # tidy data by comparisons for each guide
    df_logfc = in_dataframe.copy()
    for comp in comparisons: 
        df_logfc[comp] = in_dataframe[comp]

    return df_logfc

# calculate the negative controls (ie the mean and stdev for the non)
def calc_neg_ctrls(df_data, list_compnames, neg_ctrl_col, neg_ctrl_category): 
    """[Summary]

    This function calculates the negative control mean, stdev, upper, and lower. 
    
    Parameters
    ----------

    Returns
    ----------
    df_negctrl : pandas dataframe
        pandas dataframe of n rows and x conditions
    list_negctrlstats : list
        list of x tuples of stats
    avg_dict : dictionary
        dictionary of length x comparisons
    """

    # use negative controls to set cutoffs
    df_negctrl = df_data.loc[df_data[neg_ctrl_col].isin(neg_ctrl_category)].copy()

    list_negctrlstats = [] # list of tups of (comp, avg, sd, avg+2sd, avg-2sd)
    avg_dict = {} # dictionary of comp: avg

    for comp in list_compnames:
        # mean and stdev
        m, s = (df_negctrl[comp].mean(), df_negctrl[comp].std())
        
        # comparison, mean, 2 stdev above, 2 stdev below
        tup_comp = (comp, m, s, m + (2*s), m - (2*s))
        list_negctrlstats.append(tup_comp)
        avg_dict[comp] = m
        
    return df_negctrl, list_negctrlstats, avg_dict
