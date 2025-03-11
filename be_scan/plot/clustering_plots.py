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
import math

from pycirclize import Circos

### HELPER FUNCTIONS: PLOTTING--------------------------------------------------------------------------

def plot_PWES_heatmap(df_scaled, out_prefix, out_dir, gene_map, 
                      mask_on=True, bounds={}, 
                      sns_context='talk', cmap='RdBu_r',
                      v_min=-1, v_max=1 ):
    
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
    for gene, bound in bounds.items():
        for domain in bound: 
            if len(gene_map) <= 1: 
                # FOR ONE GENE
                # TRY START BOUND
                start = domain['start']
                while len(np.where(df_scaled.index == start)[0]) == 0 and start < df_scaled.index.max():
                    start += 1
                start = np.where(df_scaled.index == start)[0]
                # TRY END BOUND
                end = domain['end']
                while len(np.where(df_scaled.index == end)[0]) == 0 and end > df_scaled.index.min():
                    end -= 1
                end = np.where(df_scaled.index == end)[0]

            else: 
                # FOR MULTIPLE GENE
                # TRY START BOUND
                start = gene_map[gene] + str(domain['start']).zfill(4)
                while len(np.where(df_scaled.index == start)[0]) == 0 and start < df_scaled.index.max():
                    start = gene_map[gene] + str(int(start[1:])+1).zfill(4)
                start = np.where(df_scaled.index == start)[0]
                # TRY END BOUND
                end = gene_map[gene] + str(domain['end']).zfill(4)
                while len(np.where(df_scaled.index == end)[0]) == 0 and end > df_scaled.index.min():
                    end = gene_map[gene] + str(int(end[1:])-1).zfill(4)
                end = np.where(df_scaled.index == end)[0]

            if len(start) != 1 or len(end) != 1 or start[0] >= end[0]: 
                domain_name = domain['name']
                print(f'Error with domain {domain_name} in gene {gene}')
                continue

            ax.axhspan(start[0], end[0], color=domain['color'], alpha=1/8)
            ax.axvspan(start[0], end[0], color=domain['color'], alpha=1/8)

            # ADD DOMAIN NAME LABEL
            mid_x = (start[0] + end[0]) / 2
            mid_y = (start[0] + end[0]) / 2
            x_lim = ax.get_xlim()
            y_lim = ax.get_ylim()
            ax.text(x_lim[1] +5, mid_y, domain['name'], fontsize=10, ha='left', va='center', color='black', rotation=0)
            ax.text(mid_x, y_lim[1] -5, domain['name'], fontsize=10, ha='center', va='bottom', color='black', rotation=90)

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    plt.tight_layout()
    plt.savefig(f'{out_dir}/{out_prefix}_PWES_heatmap.pdf', format='pdf', bbox_inches='tight')
    plt.savefig(f'{out_dir}/{out_prefix}_PWES_heatmap.png', format='png', bbox_inches='tight')
    plt.show()
    plt.close()

def plot_clus_histogram(df_clus, 
                        out_prefix, out_dir, color_map, ):
    """
    Given df_clus, generated from cluster_pws, plot histogram of
    the number of guides in each cluster.
    """

    try:    
        df_clus['cl_new']
    except KeyError:
        raise Exception('df_clus does not contain a "cl_new" column')

    clust_counts = df_clus["cl_new"].value_counts().sort_index()
    
    plt.figure(figsize=(10, 5))
    plt.bar(clust_counts.index, clust_counts, color=[color_map[i] for i in clust_counts.index])
    
    plt.xticks(clust_counts.index)
    plt.xlabel("Cluster")
    plt.ylabel("Count")
    plt.title("Number of Guides in Cluster")
    
    plt.savefig(f'{out_dir}/{out_prefix}_cluster_histogram.pdf', format='pdf')
    plt.savefig(f'{out_dir}/{out_prefix}_cluster_histogram.png', format='png')
    plt.show()
    plt.close()

def plot_scatter_clusters(df_clus, x_col, scores_col, out_prefix, out_dir, color_map, ):
    try: 
        df_clus['cl_new']
    except KeyError:
        raise Exception('df_clus does not contain a "cl_new" column')
    
    fig, ax = plt.subplots(figsize=(15, 6))
    scatter = sns.scatterplot(data=df_clus, x=x_col, y=scores_col, hue='cl_new', palette=color_map, ax=ax)
    ax.set_xlabel("Amino Acid Position")
    ax.set_ylabel(scores_col)
    ax.set_title("Scatterplot, Colored by Cluster")

    handles, labels = ax.get_legend_handles_labels()
    ax.get_legend().remove()
    fig.savefig(f'{out_dir}/{out_prefix}_scatterplot.pdf', format='pdf')
    fig.savefig(f'{out_dir}/{out_prefix}_scatterplot.png', format='png')
    plt.show()
    plt.close(fig)

    # Save legend separately
    legend_fig, legend_ax = plt.subplots(figsize=(4, len(labels) * 0.4))
    legend_ax.axis("off")
    legend_ax.legend(handles, labels, loc='center', frameon=False)
    legend_fig.savefig(f'{out_dir}/{out_prefix}_scatterplotlegend.pdf', format='pdf')
    legend_fig.savefig(f'{out_dir}/{out_prefix}_scatterplotlegend.png', format='png')
    plt.show()
    plt.close(legend_fig)
    
def plot_clustermap(df_scaled, link, df_clusters, 
                    out_prefix, out_dir, color_list, color_clusters=True, 
                    cmap='RdBu_r', sns_context='paper', v_min=-1, v_max=1, ):

    # Plotting parameters and variables
    sns.set_context(sns_context)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['font.sans-serif'] = ['Arial']
    
    # Set colors for clusters
    if color_clusters:
        # color_list = color_list * int(
        #     np.ceil(df_clusters['cl_new'].max()/len(color_list)))
        df_colors = pd.DataFrame(index=df_scaled.index, columns=['Cluster'])
        #df_colors['Cluster'] = [color_list[i-1] for i in df_clusters['Group']]
        df_colors['Cluster'] = [color_list[i] for i in df_clusters['cl_new']]
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
    plt.savefig(f'{out_dir}/{out_prefix}_cluster_heatmap.png', format='png')
    plt.show()
    plt.close()

def plot_cluster_boxplots(df_clus, x_col, 
                          out_prefix, out_dir, color_map, ): 

    plt.figure(figsize=(12, 6))
    
    sns.boxplot(x='cl_new', y=x_col, data=df_clus, 
                palette={str(k): v for k, v in color_map.items()})

    plt.xlabel('Clusters')
    plt.ylabel(f'{x_col} Values')
    plt.title('Boxplot')
    plt.grid(True, alpha=0.5)

    plt.savefig(f'{out_dir}/{out_prefix}_cluster_boxplots.pdf', format='pdf')
    plt.savefig(f'{out_dir}/{out_prefix}_cluster_boxplots.png', format='png')
    plt.show()
    plt.close()

def plot_scatter_clusters_subplots(df_clus, x_col, scores_col, out_prefix, out_dir, color_map, cols=3):
    """
    Generates separate scatter plots for each cluster in `cl_new`,
    highlighting one cluster per plot while coloring others gray.
    """

    unique_clusters = sorted(df_clus['cl_new'].unique())  # Get unique clusters
    num_clusters = len(unique_clusters)
    rows = math.ceil(num_clusters / cols)  # Determine required number of rows

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 3))  # Create subplots
    axes = axes.flatten()  # Flatten axes array for easy indexing

    for i, cluster in enumerate(unique_clusters):
        ax = axes[i]

        # Gray out all other points
        sns.scatterplot(
            data=df_clus, x=x_col, y=scores_col, 
            color="lightgray", alpha=0.5, ax=ax, legend=False, 
        )

        # Highlight only the current cluster
        sns.scatterplot(
            data=df_clus[df_clus['cl_new'] == cluster], 
            x=x_col, y=scores_col, 
            color=color_map[cluster], 
            label=f'Cluster {cluster}', ax=ax, legend=False, 
        )

        ax.set_xlabel("Amino Acid Position")
        ax.set_ylabel(scores_col)
        ax.set_title(f"Cluster {cluster}")

    # Remove empty subplots if clusters don't perfectly fill the grid
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])  

    plt.tight_layout()
    plt.draw()
    fig.savefig(f'{out_dir}/{out_prefix}_scatterplot_subplots.pdf', format='pdf', dpi=50)
    fig.savefig(f'{out_dir}/{out_prefix}_scatterplot_subplots.png', format='png', dpi=50)
    plt.show()
    plt.close(fig)
    
def plot_scatter_clusters_complex(df_clus, x_col, scores_col, out_prefix, out_dir, color_map, gene_map):
    try: 
        df_clus['cl_new']
    except KeyError:
        raise Exception('df_clus does not contain a "cl_new" column')
    
    df_clus['chain'] = df_clus['label'].str[:1]
    df_unique_genes = df_clus['chain'].unique()
    gene_map_unique_genes = list(gene_map.values())
    common_chains = set(df_unique_genes) and set(gene_map_unique_genes)
    num_chains = len(common_chains)
    
    fig, ax = plt.subplots(1, num_chains, figsize=(15 * num_chains, 6))

    handles_list, labels_list = [], []
    for i, chain in enumerate(common_chains): 
        df_chain = df_clus[df_clus['chain'] == chain]
        scatter = sns.scatterplot(data=df_chain, x=x_col, y=scores_col, 
                                  hue='cl_new', palette=color_map, ax=ax[i])
        ax[i].set_xlabel("Amino Acid Position")
        ax[i].set_ylabel(scores_col)
        ax[i].set_title(f"Scatterplot Chain {chain} Clusters")

        handles, labels = ax[i].get_legend_handles_labels()
        ax[i].get_legend().remove()
        handles_list.append(handles)
        labels_list.append(labels)

    fig.savefig(f'{out_dir}/{out_prefix}_scatterplot.pdf', format='pdf')
    fig.savefig(f'{out_dir}/{out_prefix}_scatterplot.png', format='png')
    plt.show()
    plt.close(fig)

    # Remove duplicates while preserving order
    all_handles = sum(handles_list, [])
    all_labels = sum(labels_list, [])

    seen = set()
    unique_handles_labels = [(h, l) for h, l in zip(all_handles, all_labels) if l not in seen and not seen.add(l)]
    unique_handles, unique_labels = zip(*unique_handles_labels)

    legend_fig, legend_ax = plt.subplots(figsize=(4, len(unique_labels) * 0.4))
    legend_ax.axis("off")
    legend_ax.legend(unique_handles, unique_labels, loc='center', frameon=False)
    legend_fig.savefig(f'{out_dir}/{out_prefix}_combined_legend.pdf', format='pdf')
    legend_fig.savefig(f'{out_dir}/{out_prefix}_combined_legend.png', format='png')
    plt.show()
    plt.close(legend_fig)

def plot_scatter_clusters_subplots_complex(df_clus, x_col, scores_col, out_prefix, out_dir, color_map, gene_map):
    """
    Generates separate scatter plots for each cluster in `cl_new`,
    highlighting one cluster per plot while coloring others gray.
    """

    unique_clusters = sorted(df_clus['cl_new'].unique())  # Get unique clusters
    num_clusters = len(unique_clusters)

    df_clus['chain'] = df_clus['label'].str[:1]
    df_unique_genes = df_clus['chain'].unique()
    gene_map_unique_genes = list(gene_map.values())
    common_chains = set(df_unique_genes) and set(gene_map_unique_genes)
    num_chains = len(common_chains)

    print('Numbers', num_clusters, num_chains)
    fig, axes = plt.subplots(num_clusters, num_chains, figsize=(num_chains * 5, num_clusters * 3))  # Create subplots

    for i, cluster in enumerate(unique_clusters):
        for j, chain in enumerate(common_chains): 
            df_chain = df_clus[df_clus['chain'] == chain]
            ax = axes[i, j]

            # Gray out all other points
            sns.scatterplot(
                data=df_chain, x=x_col, y=scores_col, 
                color="lightgray", alpha=0.5, ax=ax, legend=False, 
            )
            # Highlight only the current cluster
            sns.scatterplot(
                data=df_chain[df_chain['cl_new'] == cluster], 
                x=x_col, y=scores_col, 
                color=color_map[cluster], 
                label=f'Cluster {cluster}', ax=ax, legend=False, 
            )

            ax.set_xlabel("Amino Acid Position")
            ax.set_ylabel(scores_col)
            ax.set_title(f"Chain {chain} Cluster {str(i+1)}")

    # # Remove empty subplots if clusters don't perfectly fill the grid
    # for j in range(i + 1, len(axes)):
    #     fig.delaxes(axes[j])  

    plt.tight_layout()
    plt.draw()
    fig.savefig(f'{out_dir}/{out_prefix}_scatterplot_subplots.pdf', format='pdf', dpi=50)
    fig.savefig(f'{out_dir}/{out_prefix}_scatterplot_subplots.png', format='png', dpi=50)
    plt.show()
    plt.close(fig)

def plot_PWES_heatmap_clusters(df_pwes_sorted, aas_dict, out_prefix, out_dir, 
                               mask_on=True, bounds=[], 
                               sns_context='talk', cmap='RdBu_r',
                               cols=3, v_min=-1, v_max=1 ): 

    unique_clusters = sorted(aas_dict.keys())  # Get unique clusters
    num_clusters = len(unique_clusters)
    rows = math.ceil(num_clusters / cols)  # Determine required number of rows

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 8, rows * 8))  # Create subplots
    axes = axes.flatten()  # Flatten axes array for easy indexing

    sns.set_context(sns_context)
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42    
    mpl.rcParams['font.sans-serif'] = ['Arial']

    mask = np.zeros_like(df_pwes_sorted)
    if mask_on:
        # Make a mask for lower triangle
        mask[np.tril_indices_from(mask, k=-1)] = True

    for i, clus_num in enumerate(unique_clusters):
        print(f'Cluster Number: {clus_num}')
        ax = axes[i]
        aas_list = aas_dict[clus_num]
        aas_list_floats = [float(x) for x in aas_list]

        # filtered_df = df_pwes_sorted.loc[aas_list_floats, aas_list_floats]
        df_filtered = df_pwes_sorted.copy()
        df_filtered[~df_filtered.index.isin(aas_list_floats)] = 0
        df_filtered.loc[:, ~df_filtered.columns.isin(aas_list_floats)] = 0

        sns.heatmap(data=df_filtered, ax=ax, square=True, mask=mask, cmap=cmap, 
                    xticklabels=False, yticklabels=False, vmin=v_min, vmax=v_max, 
                    cbar_kws={"shrink": .70, 'location':'left'}, )

        # LOOK AT DOMAINS
        for bound in bounds:

            # TRY START BOUND
            start = bound['start']
            while len(np.where(df_pwes_sorted.index == start)[0]) == 0 and start < df_pwes_sorted.index.max():
                start += 1
            start = np.where(df_pwes_sorted.index == start)[0]
            # TRY END BOUND
            end = bound['end']
            while len(np.where(df_pwes_sorted.index == end)[0]) == 0 and end > df_pwes_sorted.index.min():
                end -= 1
            end = np.where(df_pwes_sorted.index == end)[0]

            if start[0] >= end[0]: 
                domain_name = bound['name']
                print(f'Error with domain {domain_name}')
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

        ax.set_title(f"Cluster {clus_num}")

    # Remove empty subplots if clusters don't perfectly fill the grid
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    plt.tight_layout()
    fig.savefig(f'{out_dir}/{out_prefix}_heatmap_subplots.pdf', format='pdf')
    fig.savefig(f'{out_dir}/{out_prefix}_heatmap_subplots.png', format='png')
    plt.show()
    plt.close(fig)
    return False

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
