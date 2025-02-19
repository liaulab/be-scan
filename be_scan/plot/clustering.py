"""
Author: Christine Lee, Calvin XiaoYang Hu
Adapted from: Kevin Ngan, Nick Lue
Date: 250130

{Description: }
"""
import os

import numpy as np
import pandas as pd
import pickle

from scipy.spatial.distance import pdist, squareform
import scipy.cluster as sp_cl

from biopandas.pdb import PandasPdb
from be_scan.plot.clustering_plots import *

### HELPER FUNCTIONS: DISTANCE FACTOR (e^(-d^2 / 2t^2))--------------------------------------------------------------------------

def process_pdb(pdb_file):
    '''
    Given PDB or Alphafold .pdb file, returns dataframe with X, Y, Z coordinates and amino acid number.

    pdb_file: str, path to pdb file
    '''
    ppdb = PandasPdb().read_pdb(pdb_file)
    atom_df = ppdb.df['ATOM'] #return just atoms
    
    #gather XYZ of alpha carbons
    coord_df = atom_df.loc[atom_df["atom_name"] == "CA"] 
    df_centroids = coord_df[["residue_number", "x_coord", "y_coord", "z_coord"]] 
    df_centroids.columns = ["aa_num", "x", "y", "z"]
    
    return df_centroids

def get_pairwise_dist(df_centroids, aa_int=None):
    """
    Calculate pairwise distances from centroid coordinates.
    Returns a pandas df of pairwise distances with columns/index as AA pos
    
    df_centroids: pandas dataframe containing 4 columns of
        ['aa_num', 'x', 'y', 'z'].
    aa_int: optional tuple of (aa_min, aa_max) defining the aas to calculate pdists
        the default is None, which takes the min/max of df_centroids
    """
    
    #ARGUMENT ASSERTIONS:
    # check for correct columns in df_centroids, convert aa_num to integers
    list_cols = df_centroids.columns.tolist()
    if not all(col in list_cols for col in ['aa_num', 'x', 'y', 'z']):
        raise Exception('df_centroids is missing an essential column id')
    df_centroids['aa_num'] = df_centroids['aa_num'].astype('int64')

    # Isolate desired amino acid interval
    if aa_int is None:
        # remove unresolved residues (xyz = NaN) before finding aa min/max
        print("Removing unresolved residues...")
        df_aaint = df_centroids.loc[~df_centroids.isnull().any(axis=1)].copy()
        aa_min = df_aaint['aa_num'].min()
        aa_max = df_aaint['aa_num'].max()

    else:
        aa_min = aa_int[0]
        aa_max = aa_int[1]
        print("Removing unresolved residues...")
        df_aaint = df_centroids.loc[df_centroids['aa_num'].between(aa_min, aa_max)].copy()
        df_aaint = df_aaint.loc[~df_aaint.isnull().any(axis=1)].copy()
        if df_aaint['aa_num'].min() != aa_min:
            print('Warning! User aa_min input was ' + str(aa_min))
            print('But first resolved AA was ' + str(df_aaint['aa_num'].min()))
        if df_aaint['aa_num'].max() != aa_max:
            print('Warning! User aa_max input was ' + str(aa_max))
            print('But last resolved AA was ' + str(df_aaint['aa_num'].max()))
            
    # calculate all pairwise distances in euclidean 3d space, condense to square-form
    print("Calculating pairwise distances...")

    pairwise = pdist(df_aaint[['x','y','z']], 'euclidean')
    pairwise = squareform(pairwise)
    df_pwdist = pd.DataFrame(pairwise, index=df_aaint['aa_num'], columns=df_aaint['aa_num'])
    return df_pwdist

def gauss(distance, std):
    """
    Applies Gaussian distance function.
    """
    arg = -(distance * distance) / (2 * std * std)
    dist = np.exp(arg)
    return dist

### HELPER FUNCTIONS: ENRICHMENT SCORE CALCULATIONS-----------------------------------------------------

def hill(lfc, m, theta): 
    """
    From NZL code:
    hill function to scale the log2_fc numbers into range [0,1]
    sigmoidal func prevents highly enriched sgRNAs (e.g. jackpots) from having
    disproportionate influence vs. enriched but not jackpot sgRNAs
    input is the summed lfc (of sg1 and sg2), m, and theta
    theta controls the critical point (center), m controls steepness of function
    CLUMPS used m=3, theta=2, LSD1 used m=2, theta=3
    according to allison, m=3, theta=2 didn't work, so idk
    see the Gad Getz CLUMPS PNAS paper for reference
    """
    num = lfc**m
    denom = (lfc**m) + (theta**m)
    val = num/denom
    return val

def calculate_pw_score(df_score, scores_col, tanh_a):
    """
    Calculate pairwise sums matrix for scores.
    """

    df_pws_scores = df_score.copy()

    #Calculate pairwise sums matrix
    pw_sum = df_pws_scores[scores_col].values[:, None] + df_pws_scores[scores_col].values[None, :] # pairwise sums matrix
    df_pws_sum = pd.DataFrame(index=df_pws_scores['sgRNA_ID'], 
                                 columns=df_pws_scores['sgRNA_ID'],
                                 data= pw_sum)
    
    #Calculate parameters for tanh and zscore
    upper_tri = np.where(np.triu(np.ones(df_pws_sum.shape), k=1).astype(bool), df_pws_sum, np.nan) #get upper triangle
    pws_triu = pd.DataFrame(index=df_pws_sum.index, 
                            columns=df_pws_sum.columns, 
                            data= upper_tri)
    
    flat_pws = pd.Series([y for x in pws_triu.columns for y in pws_triu[x]], name='sum_lfc').dropna() # flatten for mean + std calc.
    df_pws = np.tanh(tanh_a * (df_pws_sum - flat_pws.mean()) / flat_pws.std())
    return df_pws


### HELPER FUNCTIONS: CALCULATE PWES-----------------------------------------------------

def signed_exp(df, exponent):
    return np.sign(df) * np.abs(df)**exponent

def calculate_pwes(df_gauss, df_pws, list_aas, pws_scaling, gauss_scaling):
    """
    Calculate PWES
    """
    df_pws.index, df_pws.columns = list_aas, list_aas #replace sgrna index with amino acid pos
    df_pws = signed_exp(df_pws, pws_scaling) * signed_exp(df_gauss.loc[list_aas, list_aas].copy(), gauss_scaling)

    #Sort by aa
    df_pws_sort = df_pws.sort_index(axis=0)
    df_pws_sort = df_pws_sort.sort_index(axis=1)

    print("PWES calculated...")

    return df_pws_sort, df_pws

def cluster_pws(df_pws, df_score, df_gauss, t, x_col):
    list_aas = df_score[df_score[x_col].isin(df_gauss.index)][x_col]
    
    #Create dendrogram
    print("Starting linkage...")
    df_clus = df_score.loc[df_score[x_col].isin(list_aas)].copy().reset_index()

    link = sp_cl.hierarchy.linkage(df_pws, method='ward', metric='euclidean', optimal_ordering=True)
    
    print("Linking complete...")
    
    df_clus['cl_new'] = sp_cl.hierarchy.fcluster(link, t=t, criterion='distance')
    #Find number of clusters
    num_clus = sorted(df_clus['cl_new'].unique())[-1]
    print(f"Number of clusters: {num_clus}")

    return df_clus, link
    
def get_clus_aa(df_clus, x_col):
    """
    Given df_clus, generated from cluster_pws, print the amino acids in each cluster.
    """
    try:    
        df_clus['cl_new']
        df_clus[x_col]
    except KeyError:
        raise Exception('df_clus does not contain required columns (cl_new or aa_pos)')

    for clus in sorted(df_clus["cl_new"].unique()):
        aas = sorted(list(set(df_clus[df_clus["cl_new"] == clus][x_col])))
        print(f'Cluster {clus} amino acids: \n{aas}')

# MAIN #

def pwes_clustering(pdb_file, scores_file, x_col, scores_col, domains_list=[], 
                    pws_scaling=1, gauss_scaling=1, 
                    gauss_std = 16, dend_t = 13.9, tanh_a=1, 
                    aa_int=None, out_prefix=None, out_dir=None):
    
    """
    Main function to run 3D clustering analysis.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Get centroids and spatial factor:
    df_centroids = process_pdb(pdb_file)
    df_pwdist = get_pairwise_dist(df_centroids, aa_int=aa_int)
    df_gauss = df_pwdist.apply(lambda x: gauss(x, gauss_std))

    # Get scores and calculate pairwise sums:
    df_scores = pd.read_csv(scores_file)
    
    sgrnaID = ["sgRNA_" + num for num in map(str, list(range(df_scores.shape[0])))] #assign id numbers 
    df_scores["sgRNA_ID"] = sgrnaID

    list_aas = df_scores[df_scores[x_col].isin(df_gauss.index)][x_col]

    df_pws_score = calculate_pw_score(df_scores, scores_col, tanh_a)

    # Calculate PWES:
    df_pwes_sorted, df_pwes_unsorted = calculate_pwes(df_gauss, df_pws_score, list_aas, 
                                                      pws_scaling, gauss_scaling)

    # Plot triangular matrix
    plot_PWES_heatmap(df_pwes_sorted, out_prefix, out_dir, 
                      bounds=domains_list)

    # Cluster PWES:
    df_clus, link = cluster_pws(df_pws = df_pwes_unsorted, 
                                df_score = df_scores, 
                                df_gauss = df_gauss,
                                t = dend_t, 
                                x_col=x_col, )

    # Plot histograms and scatterplots:
    plot_clus_histogram(df_clus, out_prefix, out_dir)
    plot_scatter_clusters(df_clus, x_col, scores_col, out_prefix, out_dir)
    plot_cluster_boxplots(df_clus, scores_col, out_prefix, out_dir)

    # Print and return clusters:
    aas_dict = get_clus_aa(df_clus, x_col)
    with open(f"{out_dir}/{out_prefix}_aas_dict.pickle", "wb") as file:
        pickle.dump(aas_dict, file)
    
    plot_PWES_heatmap_clusters(df_pwes_sorted, aas_dict, out_prefix, out_dir, 
                               bounds=domains_list)

    # Plot clustergram
    plot_clustermap(df_scaled = df_pwes_unsorted, 
                    link = link, 
                    df_clusters = df_clus, 
                    out_prefix=out_prefix, out_dir=out_dir)

    return df_pwes_sorted, df_pwes_unsorted, df_clus, aas_dict


# colors = [
#         '#beffe0', '#a2d9be', '#91c1be', '#7fa8b5', '#628397',
#         '#ead4aa', '#e8b796', '#d9ac9d', '#c99fa2', '#b08b9d',
#         '#98ff6a', '#83de6d', '#7bc17f', '#6dab8f', '#578886',
#         '#ffa296', '#f58b83', '#f6757a', '#d8677a', '#b6587c',
#         '#ff65ed', '#cf53cf', '#ae45d1', '#9a3dcc', '#7c30c8',
#         '#fffd95', '#d2d29c', '#b2b2a4', '#9898ab', '#78789e',
#         '#aec2bb', '#97a8bd', '#8494c0', '#6975b2', '#4c549b',
#         '#ffc337', '#feae34', '#e39a3b', '#be8042', '#9e6a43',
#         '#ff964e', '#ff8643', '#ff6932', '#d15842', '#b14a50',
#          ]

# EZH2_domains_list = [    # start, end, color, name
#     {'start':14, 'end':38, 
#      'color':colors[10], 'name':'SBD',},
#     {'start':38, 'end':68, 
#      'color':colors[11], 'name':'EBD',},
#     {'start':68, 'end':107, 
#      'color':colors[12], 'name':'SANBAMT',},
#     {'start':107, 'end':127, 
#      'color':colors[13], 'name':'SAL',},
#     {'start':127, 'end':157, 
#      'color':colors[14], 'name':'SRM',},
#     {'start':159, 'end':250, 
#      'color':colors[15], 'name':'SANT1',},
#     {'start':257, 'end':309, 
#      'color':colors[16], 'name':'MCSS',},
#     {'start':428, 'end':476, 
#      'color':colors[15], 'name':'SANT2',},
#     {'start':503, 'end':605, 
#      'color':colors[17], 'name':'CXC',},
#     {'start':609, 'end':728, 
#      'color':colors[18], 'name':'SET',},
# ]

# # running it with domains
# df_pwes_sorted, df_pwes_unsorted, df_clus = pwes_clustering(
#     pdb_file = 'AF_Q15910.pdb', 
#     scores_file = 'conditions_q575r_pos_mean.csv', 
#     x_col = 'xpos', 
#     scores_col  = 'DMSO-PRC2i',
#     gauss_std = 16, dend_t = 10, 
#     out_prefix='q575r_mean', 
#     out_dir='q575r_mean', 
#     domains_list=EZH2_domains_list, 
# )
# # running it without domains
# df_pwes_sorted, df_pwes_unsorted, df_clus = pwes_clustering(
#     pdb_file = 'AF_Q15910.pdb', 
#     scores_file = 'conditions_q575r_pos_mean.csv', 
#     x_col = 'xpos', 
#     scores_col  = 'DMSO-PRC2i',
#     gauss_std = 16, dend_t = 10, 
#     out_prefix='q575r_mean', 
#     out_dir='q575r_mean', 
# )
