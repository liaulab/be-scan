"""
Author: Christine Lee, Calvin XiaoYang Hu
Adapted from: Kevin Ngan, Nick Lue
Date: 250130

{Description: }
"""

import os
import pickle

from be_scan.plot.clustering_plots import *
from be_scan.plot.clustering_helpers import *
# from clustering_plots import *
# from clustering_helpers import *

# MAIN #

def pwes_clustering(
    df_scores, 
    x_col, scores_col, pdb_file, 
    gene_col='', gene_map={}, 
    domains_list={}, 
    tanh_a=1, gauss_std = 16, dend_t = 14, 
    aa_int=None, out_prefix=None, out_dir=None, 
    show_plots=True, 
    ): 
    
    """
    Main function to run 3D clustering analysis.
    """
    # BASIC CHECKS #
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    assert x_col in df_scores.columns, 'Check [x_col] in df_scores'
    assert scores_col in df_scores.columns, 'Check [scores_col] in df_scores'

    # CALCULATE SPATIAL COMPONENT FROM STRUCTURE #
    df_centroids = process_pdb(pdb_file, gene_map.values())
    df_pwdist = get_pairwise_dist(df_centroids, gene_map.values(), aa_int)
    df_gauss = df_pwdist.apply(lambda x: gauss(x, gauss_std))

    # ASSIGN IDs FOR EVERY INPUT #
    sgrnaID = ["sgRNA_" + num for num in map(str, list(range(df_scores.shape[0])))]
    df_scores["sgRNA_ID"] = sgrnaID

    df_scores[x_col] = df_scores[x_col].astype(int) # MUST BE INTEGER #
    # FOR ONLY ONE SUBUNIT #
    if len(gene_map) == 0 or gene_col == 0: 
        print('No chain(s) indicated. Automatically assigning to structure ...')
        aa_in_structure = df_centroids.aa_num.tolist()
        df_scores = df_scores[df_scores[x_col].isin(aa_in_structure)]

        # ONLY USE BASE EDITING SCORES FOR RESIDUES FOUND IN COMPLEX #
        list_aas = df_scores[x_col]
        df_pws_score = calculate_pw_score(df_scores, scores_col, tanh_a)
    # FOR A COMPLEX #
    else: 
        print('Chain(s) indicated. Mapping gene(s) to chains ...')
        assert gene_col in df_scores.columns, 'Check [gene_col] is in '

        aa_in_structure = df_centroids.label.tolist()
        df_scores = df_scores[df_scores[gene_col].isin(gene_map.keys())] # FILTER OUT GENES NOT IN MAP #
        df_scores['label'] = df_scores[gene_col].map(gene_map) + df_scores[x_col].astype(str).str.zfill(4) # FORMAT GENE-POSITION
        df_scores = df_scores[df_scores['label'].isin(aa_in_structure)]

        list_aas = df_scores['label']
        df_pws_score = calculate_pw_score(df_scores, scores_col, tanh_a)
    
    del df_centroids, df_pwdist

    # CALCULATE PWES SCORE #
    df_pwes_sorted, df_pwes_unsorted = calculate_pwes(df_gauss, df_pws_score, list_aas)

    # CLUSTER PWES #
    if len(gene_map) == 0 or gene_col == 0: cluster_xcol = x_col
    else: cluster_xcol = 'label'
    df_clus, link = cluster_pws(
        df_pws = df_pwes_unsorted, df_score = df_scores, 
        list_aas=list_aas, t=dend_t, x_col=cluster_xcol)
    
    # PLOT TRIANGLE MATRIX #
    plot_PWES_heatmap(df_pwes_sorted, gene_map, out_prefix, out_dir, 
                      bounds=domains_list, show_plots=show_plots)

    # SET CONSISTENT COLORS BETWEEN ALL PLOTS #
    unique_clusters = sorted(df_clus['cl_new'].unique())
    palette = sns.color_palette('tab20', len(unique_clusters)) # 20 COLOR PALETTE #
    color_map = {cluster: color for cluster, color in zip(unique_clusters, palette)}
    df_clus['color'] = df_clus['cl_new'].map(color_map)

    # PLOT HISTOGRAMS AND BOXPLOTS #
    plot_clust_histogram(df_clus, out_prefix, out_dir, color_map, show_plots)
    plot_clust_boxplots(df_clus, scores_col, out_prefix, out_dir, color_map, show_plots)

    if len(gene_map) == 0 or gene_col == 0: 
        plot_scatter(df_clus, x_col, scores_col, out_prefix, out_dir, color_map, show_plots)
        plot_subscatter(df_clus, x_col, scores_col, out_prefix, out_dir, color_map, show_plots)
    else: 
        plot_scatter_complex(df_clus, x_col, scores_col, out_prefix, out_dir, color_map, gene_map, show_plots)
        plot_subscatter_complex(df_clus, x_col, scores_col, out_prefix, out_dir, color_map, gene_map, show_plots)
    
    # if len(gene_map) == 0 or gene_col == 0: 
    #     plot_PWES_heatmap_clust(df_pwes_sorted, aas_dict, out_prefix, out_dir, gene_map, bounds=domains_list)

    # PLOT CLUSTERGRAM #
    plot_clustermap(df_scaled=df_pwes_sorted, link=link, df_clusters=df_clus, 
                    out_prefix=out_prefix, out_dir=out_dir, color_list=color_map, show_plots=show_plots)

    # SAVE CLUSTERS #
    aas_dict = get_clus_aa(df_clus, cluster_xcol)
    with open(f"{out_dir}/{out_prefix}_aas_dict.pickle", "wb") as file:
        pickle.dump(aas_dict, file)
        
    return df_gauss, df_pwes_sorted, df_pwes_unsorted, df_clus, aas_dict
