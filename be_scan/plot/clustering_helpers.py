"""
Author: Christine Lee, Calvin XiaoYang Hu
Adapted from: Kevin Ngan, Nick Lue
Date: 250130

{Description: }
"""

import numpy as np
import pandas as pd

from scipy.spatial.distance import pdist, squareform
import scipy.cluster as sp_cl

from biopandas.pdb import PandasPdb

### HELPER FUNCTIONS: DISTANCE FACTOR (e^(-d^2 / 2t^2))--------------------------------------------------------------------------

def process_pdb(pdb_file, chains):
    '''
    Given PDB or Alphafold .pdb file, returns dataframe with X, Y, Z coordinates and amino acid number.

    pdb_file: str, path to pdb file
    '''
    ppdb = PandasPdb().read_pdb(pdb_file)
    atom_df = ppdb.df['ATOM'] # return just atoms
    if len(chains) > 0: 
        atom_df = atom_df.loc[atom_df["chain_id"].isin(chains)]
    
    # gather XYZ of alpha carbons
    coord_df = atom_df.loc[atom_df["atom_name"] == "CA"] 
    df_centroids = coord_df[["chain_id", "residue_number", "x_coord", "y_coord", "z_coord"]] 
    df_centroids.columns = ["chain", "aa_num", "x", "y", "z"]

    if len(chains) > 0: df_centroids['label'] = df_centroids['chain'] + df_centroids['aa_num'].astype(str).str.zfill(4)
    else: df_centroids['label'] =  df_centroids['aa_num'].astype(str)
    
    return df_centroids

def get_pairwise_dist(df_centroids, chains, aa_int=None):
    """
    Calculate pairwise distances from centroid coordinates.
    Returns a pandas df of pairwise distances with columns/index as AA pos
    
    df_centroids: pandas dataframe containing 4 columns of
        ['aa_num', 'x', 'y', 'z'].
    aa_int: optional tuple of (aa_min, aa_max) defining the aas to calculate pdists
        the default is None, which takes the min/max of df_centroids
    """
    
    df_centroids['aa_num'] = df_centroids['aa_num'].astype('int64')

    # Isolate desired amino acid interval
    if aa_int is None:
        # remove unresolved residues (xyz = NaN) before finding aa min/max
        print("Removing unresolved residues...")
        df_aaint = df_centroids.loc[~df_centroids.isnull().any(axis=1)].copy()

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
    
    if len(chains) > 0: df_pwdist = pd.DataFrame(pairwise, index=df_aaint['label'], columns=df_aaint['label'])
    else: df_pwdist = pd.DataFrame(pairwise, index=df_aaint['aa_num'], columns=df_aaint['aa_num'])
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

    # CALCULATE PAIRWISE SCORE OF EDITING DATA #
    pw_sum = df_pws_scores[scores_col].values[:, None] + df_pws_scores[scores_col].values[None, :] # PAIRWISE SUM MATRIX #
    df_pws_sum = pd.DataFrame(index=df_pws_scores['sgRNA_ID'], 
                              columns=df_pws_scores['sgRNA_ID'], data=pw_sum)
    
    # CALC PARAMS FOR TANH OF Z-SCORE #
    upper_tri = np.where(np.triu(np.ones(df_pws_sum.shape), k=1).astype(bool), 
                         df_pws_sum, np.nan) # GET UPPER TRIANGLE #
    pws_triu = pd.DataFrame(index=df_pws_sum.index, 
                            columns=df_pws_sum.columns, data=upper_tri)
    
    flat_pws = pd.Series([y for x in pws_triu.columns for y in pws_triu[x]], name='sum_lfc').dropna() # FLATTEN FOR MEAN + STD CALC #
    df_pws = np.tanh(tanh_a * (df_pws_sum - flat_pws.mean()) / flat_pws.std()) # SCALED TANH FUNCTION #
    # COULD ALSO USE HILL FUNCTION HERE #
    return df_pws

### HELPER FUNCTIONS: CALCULATE PWES-----------------------------------------------------

def calculate_pwes(df_gauss, df_pws, list_aas):
    """
    Calculate PWES
    """
    df_pws.index, df_pws.columns = list_aas, list_aas
    # MAKE SURE THE INDEX OF df_pws AND df_gauss ARE THE SAME #
    df_pws = df_pws * df_gauss.loc[list_aas, list_aas].copy()

    # SORT BY AA #
    df_pws_sort = df_pws.sort_index(axis=0)
    df_pws_sort = df_pws_sort.sort_index(axis=1)
    print("PWES calculated ...")

    return df_pws_sort, df_pws

def cluster_pws(df_pws, df_score, list_aas, t, x_col):
    
    # Create dendrogram
    print("Starting linkage ...")
    df_clus = df_score.loc[df_score[x_col].isin(list_aas)].copy().reset_index()
    # df_clus['label_encoded'] = df_clus['label'].astype('category').cat.codes ###
    link = sp_cl.hierarchy.linkage(df_pws, method='ward', metric='euclidean', optimal_ordering=True)
    print("Linking complete ...")
    
    df_clus['cl_new'] = sp_cl.hierarchy.fcluster(link, t=t, criterion='distance')
    # Find number of clusters
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

    aas_dict = {}
    for clus in sorted(df_clus["cl_new"].unique()):
        df_subset = df_clus[df_clus["cl_new"] == clus][x_col]
        aas = sorted([x for x in list(set(df_subset))])
        print(f'Cluster {clus} amino acids: \n{aas}')
        aas_dict[clus] = aas
    return aas_dict

def shuffle_pwes(df_gauss, df_pws, list_aas, nrand=1000):
    """
    Randomize PWES Score
    """
    df_pws.index, df_pws.columns = list_aas, list_aas
    df_pws_list = []
    np.random.seed(0)

    for i in range(nrand): 
        np.random.seed(i)
        permutation = np.random.permutation(df_pws.shape[0])
        df_shuffled_sample = df_pws.iloc[permutation, permutation] # SHUFFLE ROWS AND COLS #
        df_shuffled_sample.index, df_shuffled_sample.columns = list_aas, list_aas

        df_shuffled_sample *= df_gauss.loc[list_aas, list_aas].copy()
        df_pws_list.append(df_shuffled_sample)

    # COLLAPSE RANDOMIZATIONS INTO ONE MATRIX #
    result_df = sum(df_pws_list) / len(df_pws_list)
    # SORT BY AA #
    df_pws_sort = result_df.sort_index(axis=0).sort_index(axis=1)

    print("Randomized PWES calculated.")
    return df_pws_sort, result_df
