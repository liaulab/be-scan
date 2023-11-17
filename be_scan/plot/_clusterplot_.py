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

#%% plot_clustermap() - plots the hierarchical clustered heatmap and dendrogram
# Made by KCN, modified by NZL

def plot_clustermap(df_scaled, link, df_clusters, num_clusters, out_prefix,
                    cmap='RdBu_r', sns_context='talk', v_min=-1, v_max=1,
                    color_clusters=True,
                    color_list=list(sns.color_palette('deep').as_hex())):

    # Plotting parameters and variables
    sns.set_context(sns_context)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Set colors for clusters
    if color_clusters:
        color_list = color_list * int(
            np.ceil(df_clusters['Group'].max()/len(color_list)))
        df_colors = pd.DataFrame(index=df_scaled.index, columns=['Cluster'])
        df_colors['Cluster'] = [color_list[i-1] for i in df_clusters['Group']]
        rcolor = df_colors['Cluster'].copy()
    else:
        rcolor = None

    fig = sns.clustermap(df_scaled, row_linkage=link, col_linkage=link,
                         cmap=cmap, vmin=v_min, vmax=v_max, row_colors=rcolor, 
                         xticklabels=False, yticklabels=False, center=0)
    
    # Check clusters ordered correctly. Returned list should be monotonic.
    cluster_check = [df_clusters['Group'][i] for i in fig.dendrogram_row.reordered_ind]
    print(all(cluster_check[i]<=cluster_check[i+1] for i in range(len(cluster_check)-1)))
    
    fig.ax_col_dendrogram.set_visible(False)
    ax = fig.ax_heatmap
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_aspect('equal')
    # plt.title("Clustered Heatmap")
    clusterlist = df_clusters['Group'].value_counts().sort_index()
    temp = 0
    for i in clusterlist.iloc[:-1]:
        temp = temp + i
        ax.axhline(y=temp, color='k', lw=1, ls='--')
    for edge,spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
    
    plt.savefig(out_prefix + '_' + str(num_clusters) + 'cluster_heatmap.pdf',
                format='pdf')
    plt.close()
