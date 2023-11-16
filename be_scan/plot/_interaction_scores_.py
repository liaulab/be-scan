"""
Author: Calvin XiaoYang Hu
Adapted from: Nicholas Lue - NZL10196-3D-Clustering-v7.py Created on Mon Jun  1 19:03:10 2020
Date: 231116

{Description: some base pair to amino acid translation functions}
"""

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF

def plot_interaction_scores(df_plot, score_col, enr_col, marker_color,
                            list_negctrlstats, plot_name, out_prefix,
                            hue_col='Mut_type', hue_order=None, palette=None,
                            axis_title=None,
                            sns_context='paper', sns_palette='deep'):
    '''
    Generate plots summarizing interaction scores for each guide.
    Left: scatterplot showing enrichments vs. interaction scores.
    Right: cumulative distribution function of interaction scores.
    '''
    
    # Set seaborn/mpl parameters
    sns.set_context(sns_context)
    sns.set_palette(sns_palette)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Initialize plot
    fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(8,4))

    # Left plot: scatterplot of summed PWES vs. lfc val, each point is a guide
    sns.scatterplot(x=enr_col, y=score_col, data=df_plot, ax=ax1,
                    alpha=0.8, hue=hue_col, hue_order=hue_order,
                    palette=palette, legend=False,
                    linewidth=1, edgecolor='black', s=50)
    # Plot vertical dashed lines denoting avg +/- 2 stdev of controls
    tup_plot = [tup for tup in list_negctrlstats if tup[0] == enr_col][0]
    ax1.axvline(x=2*tup_plot[2], color='k', ls='--', lw=1)
    ax1.axvline(x=-2*tup_plot[2], color='k', ls='--', lw=1)
    # Aesthetic adjustments
    for edge,spine in ax1.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
    ax1.set_xlabel('d9-pos normalized log2 enrichment')
    if axis_title == None:
        ax1.set_ylabel('Interaction score')
    else:
        ax1.set_ylabel(axis_title)
    
    # Right plot: cumulative distribution plot of summed PWES for each guide
    pwes_ecdf = ECDF(df_plot[score_col])
    sns.lineplot(x=pwes_ecdf.x, y=pwes_ecdf.y, ax=ax2, color=marker_color,
                 linewidth=2)
    # Aesthetic adjustments
    for edge,spine in ax2.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
    if axis_title == None:
        ax2.set_xlabel('Interaction Score')
    else:
        ax2.set_xlabel(axis_title)
    ax2.set_ylabel('Cumulative Probability')
    
    # Adjust dimensions
    plt.tight_layout()
    
    # Save figure and close
    file_name = '_'.join([out_prefix, plot_name, 'interaction_scores.pdf'])
    plt.savefig(file_name, format='pdf', bbox_inches='tight')
    plt.close()
