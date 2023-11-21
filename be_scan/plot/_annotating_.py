"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196_Screen_Analysis_v9b.py Created on Fri May 29 03:00:39 2020
Date: 231116

{Description: }
"""

# ColorBrewer2, 9 data classes, qualitative, 4th color scheme hex codes
color_list = ['#fb8072', '#80b1d3', '#fdb462', '#b3de69', '#fccde5', 
            '#d9d9d9', '#8dd3c7', '#ffffb3', '#bebada', '#bc80bd']
# Define lists to use for setting plotting parameters
list_muttypes = ['Nonsense', 'Missense', 'Silent', 'Non-exon', 'Splice', 
                'No_C/Exon', 'No_C/Non-exon', 'Control']

def norm_to_intergenic_ctrls(in_dataframe, comparisons, avg_dict, y_column): 
    # normalize data to intergenic controls
    # perform normalization
    for comp in comparisons:
        in_dataframe[comp] = in_dataframe[comp].sub(avg_dict[comp])
    # tidy data by comparisons for each guide
    df_logfc = in_dataframe.copy()
    for comp in comparisons: 
        df_logfc[comp+'_'+y_column] = in_dataframe[comp]
    return df_logfc

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



# # Add extra label for control guides to subtype as intergenic, non-targeting, essential. 
# def annotate_submuttype(mutType, geneID):
#     if mutType == 'Control':
#         try: 
#             return control_subtypes[geneID]
#         except KeyError: 
#             print('Invalid type for control guide')
#     else: 
#         return mutType

# # For plotting domain vs. not domain
# def annotate_in_domain(domain, domains_list):
#     return 'in domain' if domain in domains_list else 'not in domain'
