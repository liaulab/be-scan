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

import scipy as sp
import scipy.spatial.distance as dist
import scipy.cluster as sp_cl

# from biopandas.pdb import PandasPdb

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

    pairwise = dist.pdist(df_aaint[['x','y','z']], 'euclidean')
    pairwise = dist.squareform(pairwise)
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


def calculate_pw_score(df_score, scores_col):
    """
    Calculate pairwise sums matrix for scores
    """

    df_pws_scores = df_score.copy()

    #Calculate pairwise sums matrix
    pw_sum = df_pws_scores[scores_col].values[:, None] + df_pws_scores[scores_col].values[None, :] # pairwise sums matrix
    df_pws_sum = pd.DataFrame(index=df_pws_scores['sgRNA_ID'], 
                                 columns=df_pws_scores['sgRNA_ID'],
                                 data= pw_sum)
    
    print(df_pws_sum)
    
    #Calculate parameters for tanh and zscore
    upper_tri = np.where(np.triu(np.ones(df_pws_sum.shape), k=1).astype(bool), df_pws_sum, np.nan) #get upper triangle
    pws_triu = pd.DataFrame(index=df_pws_sum.index, 
                            columns=df_pws_sum.columns, 
                            data= upper_tri)
    
    print(pws_triu)
    flat_pws = pd.Series([y for x in pws_triu.columns for y in pws_triu[x]], name='sum_lfc').dropna() # flatten for mean + std calc.
    print(flat_pws)
    df_pws = np.tanh((df_pws_sum - flat_pws.mean()) / flat_pws.std())
    print(df_pws)
    return df_pws


### HELPER FUNCTIONS: CALCULATE PWES-----------------------------------------------------

def calculate_pwes(df_gauss, df_pws, list_aas):
    """
    Calculate PWES
    """
    print(df_pws)
    print(df_gauss.loc[list_aas, list_aas])
    df_pws.index, df_pws.columns = list_aas, list_aas #replace sgrna index with amino acid pos 
    print(df_pws)
    df_pws = df_pws * df_gauss.loc[list_aas, list_aas].copy()

    #Sort by aa
    df_pws_sort = df_pws.sort_index(axis=0)
    df_pws_sort = df_pws_sort.sort_index(axis=1)

    print(df_pws_sort)

    print("PWES calculated...")

    return df_pws_sort

def cluster_pws(df_pws, df_score, df_gauss, scores_col, t, out_prefix = None):
    list_aas = df_score[df_score['aa_pos'].isin(df_gauss.index)]['aa_pos']
    
    print("Starting linkage...")
    df_clus = df_score.loc[df_score['aa_pos'].isin(list_aas)].copy().reset_index()
    link = sp_cl.hierarchy.linkage(df_pws, method='ward', metric='euclidean', optimal_ordering=True)
    
    print("Linking complete...")
    
    df_clus['cl_new'] = sp_cl.hierarchy.fcluster(link, t=t, criterion='distance')
    df_sorted = df_clus.sort_values(by='cl_new')
    
    num_clus = sorted(df_clus['cl_new'].unique())[-1]
    print(f"Number of clusters: {num_clus}")
    
    clust_counts = df_clus["cl_new"].value_counts().sort_index()
    cluster_histogram = plt.figure()
    plt.bar(clust_counts.index, clust_counts)
    plt.xticks(clust_counts.index)
    plt.xlabel("Cluster")
    plt.ylabel("Count")
    plt.title("Number of Guides in Cluster")
    plt.savefig(out_prefix + 'cluster_histogram.pdf', format='pdf')
    plt.show()
    
    scatter_clusters = plt.figure(figsize = (15, 9))
    sns.scatterplot(data = df_clus, x = 'aa_pos', y = scores_col, hue = 'cl_new', palette = 'bright')
    plt.xlabel("Amino Acid Position")
    plt.ylabel(scores_col)
    plt.title("Scatterplot, Colored by Cluster")
    plt.savefig(out_prefix + 'scatterplot.pdf', format='pdf')
    plt.show()
    
    plt.close()
    print("Scatterplot complete...")
    for clus in sorted(df_clus["cl_new"].unique()):
        aas = sorted(list(set(df_clus[df_clus["cl_new"] == clus]["aa_pos"])))
        print(f'Cluster {clus} amino acids: \n{aas}')
    
    plot_clustermap(df_scaled = df_pws, link = link, df_clusters = df_clus, num_clusters = num_clus, out_prefix = out_prefix)


### HELPER FUNCTIONS: PLOTTING--------------------------------------------------------------------------
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
    plt.savefig(out_prefix+'PWES_heatmap.pdf', format='pdf', bbox_inches='tight')
    plt.show()
    plt.close()
    

def plot_clustermap(df_scaled, link, df_clusters, num_clusters, out_prefix,
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
    plt.savefig(out_prefix + 'cluster_heatmap.pdf', format='pdf')
    plt.show()
    
    plt.close()


def pwes_clustering(
        pdb_file, scores_file, scores_col,
        norm_type = "tanh",
        gauss_std = 16, dend_t = 13.9, 
        aa_int=None, out_prefix=None):
    
    """
    Main function to run 3D clustering analysis.
    """
    if not os.path.exists(out_prefix):
        os.makedirs(out_prefix)

    #Get centroids and spatial factor:
    df_centroids = process_pdb(pdb_file)
    df_pwdist = get_pairwise_dist(df_centroids, aa_int=aa_int)
    df_gauss = df_pwdist.apply(lambda x: gauss(x, gauss_std))

    #Get scores and calculate pairwise sums:
    df_scores = pd.read_csv(scores_file)
    
    print(df_scores)
    sgrnaID = ["sgRNA_" + num for num in map(str, list(range(df_scores.shape[0])))] #assign id numbers 
    df_scores["sgRNA_ID"] = sgrnaID

    list_aas = df_scores[df_scores['aa_pos'].isin(df_gauss.index)]['aa_pos']
    print(df_scores["aa_pos"])
    print(df_scores)
    df_pws_score = calculate_pw_score(df_scores, scores_col)

    #Calculate PWES:
    #if(norm_type == "tanh"):
    df_pwes = calculate_pwes(df_gauss, df_pws_score, list_aas)

    #Plot triangular matrix
    plot_PWES_heatmap(df_pwes, 
                      out_prefix)

    #Cluster PWES:
    cluster_pws(df_pws = df_pwes, 
            df_score = df_scores, 
            df_gauss = df_gauss,
            scores_col = scores_col,
            t = 13,
            out_prefix = out_prefix)

#TEST
# main(pdb_file = "./data/Alphafold/AF_Q15022.pdb",
#      scores_file = "./data/clean_data/SUZ12_clean_data_WT_ABE.csv",
#      scores_col = "PRC2i_minus_DMSO_Z",
#      out_prefix = "./outputs/test_heatmaps/")
