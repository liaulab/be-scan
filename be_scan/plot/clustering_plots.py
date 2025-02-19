"""
Author: Christine Lee, Calvin XiaoYang Hu
Adapted from: Kevin Ngan, Nick Lue
Date: 250130

{Description: }
"""

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from pycirclize import Circos

### HELPER FUNCTIONS: PLOTTING--------------------------------------------------------------------------

def plot_clus_histogram(df_clus, 
                        out_prefix, out_dir, ):
    """
    Given df_clus, generated from cluster_pws, plot histogram of
    the number of guides in each cluster.
    """

    try:    
        df_clus['cl_new']
    except KeyError:
        raise Exception('df_clus does not contain a "cl_new" column')

    
    clust_counts = df_clus["cl_new"].value_counts().sort_index()
    plt.bar(clust_counts.index, clust_counts)
    plt.xticks(clust_counts.index)
    plt.xlabel("Cluster")
    plt.ylabel("Count")
    plt.title("Number of Guides in Cluster")
    plt.legend(loc="upper left", bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    plt.savefig(f'{out_dir}/{out_prefix}_cluster_histogram.pdf', format='pdf')
    plt.close()
    plt.show()

def plot_scatter_clusters(df_clus, x_col, scores_col, 
                          out_prefix, out_dir, ):
    try: 
        df_clus['cl_new']
    except KeyError:
        raise Exception('df_clus does not contain a "cl_new" column')
    plt.figure(figsize=(15, 6))
    sns.scatterplot(data = df_clus, x = x_col, y = scores_col, hue = 'cl_new', palette = 'bright')
    plt.xlabel("Amino Acid Position")
    plt.ylabel(scores_col)
    plt.title("Scatterplot, Colored by Cluster")
    plt.savefig(f'{out_dir}/{out_prefix}_scatterplot.pdf', format='pdf')
    plt.show()
    plt.close()

def plot_PWES_heatmap(df_scaled, out_prefix, out_dir, 
                      mask_on=True, bounds=[],
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
    fig, ax = plt.subplots(figsize=(10, 8))
        
    sns.heatmap(data=df_scaled, ax=ax, square=True, mask=mask, cmap=cmap,
                xticklabels=False, yticklabels=False, vmin=v_min, vmax=v_max,
                cbar_kws={"shrink": .70, 'location':'left'}, )
    
    for edge, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
    
    # LOOK AT DOMAINS
    for bound in bounds:

        # TRY START BOUND
        start = bound['start']
        while len(np.where(df_scaled.index == start)[0]) == 0 and start < df_scaled.index.max():
            start += 1
        start = np.where(df_scaled.index == start)[0]
        # TRY END BOUND
        end = bound['end']
        while len(np.where(df_scaled.index == end)[0]) == 0 and end > df_scaled.index.min():
            end -= 1
        end = np.where(df_scaled.index == end)[0]

        if start[0] >= end[0]: 
            name_temp = bound['name']
            print(f'Error with domain {name_temp}')
            continue

        ax.axhspan(start[0], end[0], color=bound['color'], alpha=1/8)
        ax.axvspan(start[0], end[0], color=bound['color'], alpha=1/8)

        # ADD DOMAIN NAME LABEL
        mid_x = (start[0] + end[0]) / 2
        mid_y = (start[0] + end[0]) / 2
        x_lim = ax.get_xlim()
        y_lim = ax.get_ylim()
        ax.text(x_lim[1] +5, mid_y, bound['name'], fontsize=10, ha='left', va='center', color='black', rotation=0)
        ax.text(mid_x, y_lim[1] -5, bound['name'], fontsize=10, ha='center', va='bottom', color='black', rotation=90)

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.tight_layout()
    plt.savefig(f'{out_dir}/{out_prefix}_PWES_heatmap.pdf', format='pdf', bbox_inches='tight')
    plt.show()
    plt.close()
    
def plot_clustermap(df_scaled, link, df_clusters, 
                    out_prefix, out_dir, 
                    cmap='RdBu_r', sns_context='paper', v_min=-1, v_max=1,
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
            np.ceil(df_clusters['cl_new'].max()/len(color_list)))
        df_colors = pd.DataFrame(index=df_scaled.index, columns=['Cluster'])
        #df_colors['Cluster'] = [color_list[i-1] for i in df_clusters['Group']]
        df_colors['Cluster'] = [color_list[i-1] for i in df_clusters['cl_new']]
        rcolor = df_colors['Cluster'].copy()
    else:
        rcolor = None

    fig = sns.clustermap(df_scaled, row_linkage=link, col_linkage=link,
                         cmap=cmap, vmin=v_min, vmax=v_max, row_colors=rcolor, 
                         xticklabels=False, yticklabels=False, center=0, cbar_pos=(1.01, 0.2, 0.03, 0.5))
    
    # Check clusters ordered correctly. Returned list should be monotonic.
    cluster_check = [df_clusters['cl_new'][i] for i in fig.dendrogram_row.reordered_ind]
    print(all(cluster_check[i]<=cluster_check[i+1] for i in range(len(cluster_check)-1)))
    
    fig.ax_col_dendrogram.set_visible(False)
    ax = fig.ax_heatmap
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_aspect('equal')
    #plt.title("Clustered Heatmap")
    clusterlist = df_clusters['cl_new'].value_counts().sort_index()
    temp = 0
    for i in clusterlist.iloc[:-1]:
        temp = temp + i
        ax.axhline(y=temp, color='k', lw=1, ls='--')
    for edge,spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_color('k')
    plt.savefig(f'{out_dir}/{out_prefix}_cluster_heatmap.pdf', format='pdf')
    plt.show()
    plt.close()

def plot_cluster_boxplots(df_clus, x_col, 
                          out_prefix, out_dir, 
                          ): 

    # Set figure size
    plt.figure(figsize=(12, 6))

    # Create the boxplot
    sns.boxplot(x='cl_new', y=x_col, data=df_clus)

    # Customize plot
    plt.xlabel('Clusters')
    plt.ylabel(f'{x_col} Values')
    plt.title('Boxplot')
    plt.grid(True, alpha=0.5)

    # Show plot
    plt.savefig(f'{out_dir}/{out_prefix}_cluster_boxplots.pdf', format='pdf')
    plt.show()
    plt.close()

from pycirclize import Circos

def plot_pwes_circos(
        df, out_prefix, out_dir, 
        ignore_range=0, width=30, 
        EZH2_domains_list=[]): 

    aa_pos_list = df.columns
    sectors = dict(zip(aa_pos_list, [width]*len(aa_pos_list))) # SET WIDTH SO THERES SPACE FOR MORE CONNECTIONS
    circos = Circos(sectors)
    EZH2_domains_bools = [True]*len(EZH2_domains_list)
    # EZH2_domains_dict = dict(zip(EZH2_domains_list, [True]*len(EZH2_domains_list)))

    # SET UP THE COLORS ACCORDING TO DOMAIN ANNOT DICTIONARY #
    for i, (sector, res_num) in enumerate(zip(circos.sectors, aa_pos_list)): 
        track = sector.add_track((100, 104))

        # IF NO ANNOT, COLOR WHITE AND ANNOT NUM EVERY 20 RESIDUES
        if len(EZH2_domains_list) == 0: 
            track.axis(fc='white')
            if i % 20 == 0: 
                sector.text(res_num, r=110, size=10)
            continue
        
        # IF THERE IS DOMAIN ANNOT
        for i, (item, bool) in enumerate(zip(EZH2_domains_list, EZH2_domains_bools)): 
            if item['start'] <= int(res_num) <= item['end']: # IF RESIDUE IS WITHIN DOMAIN
                track.axis(fc=item['color']) # RECOLOR

                # ANNOT DOMAIN NUMBERS
                if int(res_num) == item['start'] or int(res_num) == item['end']: 
                    sector.text(res_num, r=110, size=10) # ADD NUM ANNOT

                midpoint = int((item['start'] + item['end'])/2)
                # ANNOT DOMAIN NAME
                if EZH2_domains_bools[i] and midpoint-2 < int(res_num) < midpoint+2: 
                    ### this likely misses some name annotations if there i no residue with a +-2 range
                    sector.text(item['name'], r=110, size=15) # ADD LABEL ANNOT
                    EZH2_domains_bools[i] = False # SKIP THIS DOMAIN IN THE FUTURE

                break # IF RES WITHIN DOMAIN, WE CAN SKIP REST OF DOMAINS TO NEXT RES

    matrix = df.values.tolist()
    matrix_abs = np.abs(matrix)
    matrix_max = np.max(matrix_abs)
    matrix_alpha = matrix_abs / matrix_max

    int1, int2 = 5, 6
    for i in range(len(aa_pos_list)): 
        for j in range(len(aa_pos_list)): 
            if i <= j or j - ignore_range < i < j + ignore_range:
                # ONLY HALF OF THE MATRIX, AND IGNORE CLOSE RANGE INTERACTIONS
                continue

            val = matrix[i][j]
            # SKIP ZERO VALS #
            if matrix[i][j] == 0.0: 
                ### could also change this to a cmap colorbar instead of binary
                continue

            if val < 0.0: color = 'royalblue' # blue
            if val > 0.0: color = 'orangered' # red
            alpha_ratio_dec = matrix_alpha[i][j]

            # MANUALLY SETTING LINK WITH START, END, COLOR, TRANSPARENCY
            circos.link((aa_pos_list[i], int1, int2), (aa_pos_list[j], int1, int2), 
                        color=color, alpha=alpha_ratio_dec*.3)
            # MANUALLY SET HOW TO DRAW CIRCOS PLOT
            int1, int2 = (int1 + 1) % width, (int2 + 1) % width

    fig = circos.plotfig(dpi=300, figsize = (15, 15))
    fig.savefig(f'{out_dir}/{out_prefix}_pwes_circos.pdf', format='pdf')
