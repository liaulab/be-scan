
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
