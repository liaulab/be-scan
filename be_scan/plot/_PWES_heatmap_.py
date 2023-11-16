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

#%% plot_pw_heatmap() - plots the scaled LFCs as a Hi-C type heatmap

# Make heatmap for PWES scores, guides still arranged by Clust_site.
# Based on KCN function, modified by NZL
def plot_PWES_heatmap(df_scaled, out_prefix, mask_on=True, bounds=[],
                      mark_bounds=True, sns_context='talk', cmap='RdBu_r',
                      v_min=-1, v_max=1):
    
    # Plotting parameters and variables
    sns.set_context(sns_context)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42    
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Make mask; all values set to False initially (mask off)
    mask = np.zeros_like(df_scaled)
    if mask_on:
        # Make a mask for lower triangle
        mask[np.tril_indices_from(mask, k=-1)] = True
    
    # Generate heatmap
    fig, ax = plt.subplots(figsize=(5,4))
    sns.heatmap(data=df_scaled, ax=ax, square=True, mask=mask, cmap=cmap,
                xticklabels=False, yticklabels=False, vmin=v_min, vmax=v_max,
                cbar_kws={"shrink": .70})
    for edge,spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
    if mark_bounds:
        for i in bounds:
            temp = np.where(df_scaled.index == i)[0]
            ax.axhline(y=temp[0], color='k', ls='--')
            ax.axvline(x=temp[0], color='k', ls='--')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.tight_layout()
    plt.savefig(out_prefix+'_PWES_heatmap.pdf', format='pdf', bbox_inches='tight')
    plt.close()
