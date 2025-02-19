"""
Author: Christine Lee, Calvin XiaoYang Hu
Adapted from: Kevin Ngan, Nick Lue
Date: 250130

{Description: }
"""
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.spatial.distance import pdist, squareform
import scipy.cluster as sp_cl

from biopandas.pdb import PandasPdb

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

def pwes_clustering(pdb_file, scores_file, x_col, scores_col, domains_list=[], 
                    norm_type = "tanh", pws_scaling=1, gauss_scaling=1, 
                    gauss_std = 16, dend_t = 13.9, tanh_a=1, 
                    aa_int=None, out_prefix=None, out_dir=None):
    
    """
    Main function to run 3D clustering analysis.
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    #Get centroids and spatial factor:
    df_centroids = process_pdb(pdb_file)
    df_pwdist = get_pairwise_dist(df_centroids, aa_int=aa_int)
    df_gauss = df_pwdist.apply(lambda x: gauss(x, gauss_std))

    #Get scores and calculate pairwise sums:
    df_scores = pd.read_csv(scores_file)
    
    sgrnaID = ["sgRNA_" + num for num in map(str, list(range(df_scores.shape[0])))] #assign id numbers 
    df_scores["sgRNA_ID"] = sgrnaID

    list_aas = df_scores[df_scores[x_col].isin(df_gauss.index)][x_col]

    df_pws_score = calculate_pw_score(df_scores, scores_col, tanh_a)

    #Calculate PWES:
    df_pwes_sorted, df_pwes_unsorted = calculate_pwes(df_gauss, df_pws_score, list_aas, 
                                                      pws_scaling, gauss_scaling)

    #Plot triangular matrix
    plot_PWES_heatmap(df_pwes_sorted, out_prefix, out_dir, 
                      bounds=domains_list)

    #Cluster PWES:
    df_clus, link = cluster_pws(df_pws = df_pwes_unsorted, 
                                df_score = df_scores, 
                                df_gauss = df_gauss,
                                t = dend_t, 
                                x_col=x_col, )

    #Plot histograms and scatterplots:
    plot_clus_histogram(df_clus, out_prefix, out_dir)
    plot_scatter_clusters(df_clus, x_col, scores_col, out_prefix, out_dir)
    plot_cluster_boxplots(df_clus, scores_col, out_prefix, out_dir)

    #Print clusters:
    get_clus_aa(df_clus, x_col)

    #Plot clustergram
    plot_clustermap(df_scaled = df_pwes_unsorted, 
                    link = link, 
                    df_clusters = df_clus, 
                    out_prefix=out_prefix, out_dir=out_dir)

    return df_pwes_sorted, df_pwes_unsorted, df_clus


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
