"""
Author: Calvin XiaoYang Hu
Adapted from: Kevin Ngan - DNMT_Figure2_analysis_functions.py 
Date: 240516

{Description: Helper functions for clustering.py
              get_centroids calculates the distances between all combination of amino acids
              get_pairwise_dist calculates the 
              }
"""

import numpy as np
import pandas as pd
import scipy.spatial.distance as sp_sd

from pymol import cmd
from pymol import stored


def get_centroids(pdb_id, 
                  aa_sel=None, list_aas=None, save_centroids=True,
                  out_csv=None, return_df=True
                  ):
    
    """
    This function gets all centroids from a pymol structure and outputs as csv

    Parameters
    ----------
    pdb_id : str
        pdb id of the desired structure
    aa_sel : str
        str in pymol syntax for the desired aa selection (NO PARENTHESES)
        The default is pdb_id + chain A + polymer.protein
    list_aas : list
        list of all amino acids to find centroids for.
        By default uses first to last AA in aa_sel w/ resolved alpha carbons
    save_centroids : TYPE, optional
        choose whether to save the centroid df as csv. The default is True.
    out_csv : str
        name of the output csv file. The default is pdb_id + _centroids.csv
    return_df : bool
        choose whether or not to return the centroid df (default is True)
    """

    # fetch pdb structure
    cmd.fetch(pdb_id)
    # if pdb_id == '4wxx':
    #    cmd.remove('chain B')
    # default selection for centroids is the entire crystal structure
    # if a custom selection is defined, MUST be in pymol syntax
    if aa_sel is None:
        aa_sel = '%s and chain A and polymer.protein' % pdb_id
    if list_aas is None:
        # make list to hold all amino acids in the aa_sel
        stored.aa=[]
        # iterate thru AA (alpha Cs), append (resi num, resi name) to stored.list
        # pymol syntax iterate (selection), expression
        cmd.iterate('(' + aa_sel + ' and name ca)', 'stored.aa.append((int(resi),resn))')
        list_aas = [x for x in range(min(stored.aa)[0], max(stored.aa)[0] + 1)]

    # iterate through residues and add centroids to df as (x,y,z)
    # determine centroid by average of atoms(x,y,z) coords
    centroid_list = []
    for aa in list_aas:
        coords = cmd.get_coords('(' + aa_sel + ' and resi ' + str(aa) + ')')
        # if residue is unresolved, get_coords returns NoneType
        # therefore, skip unresolved residues (fill w/ None)
        if coords is None:
            centroid = (aa, None, None, None)
            centroid_list.append(centroid)
        else:
            xyz = np.mean(coords, axis=0)
            centroid = (aa, xyz[0], xyz[1], xyz[2])
            centroid_list.append(centroid)

    # turn list into pandas dataframe of xyz centroid coords
    df_centroids = pd.DataFrame(data=centroid_list,
                                columns=['aa_num', 'x', 'y', 'z'])
    # pymol sometimes makes numbers into strings -- make sure they are int
    df_centroids['aa_num'] = df_centroids['aa_num'].astype(int)
    num_resolved = df_centroids.loc[~df_centroids['x'].isna()].shape[0]

    if out_csv is None:
        out_csv = pdb_id + '_centroids.csv'
    # save centroid coordinates as csv
    if save_centroids:
        df_centroids.to_csv(out_csv, index=False)
    # print statistics
    print('# residues in list_aas: ' + str(len(list_aas)))
    print('# centroids in df_centroids: ' + str(df_centroids.shape[0]))
    print('# centroids resolved: ' + str(num_resolved))
    if return_df:
        print('get_centroids done for ' + pdb_id)
        return df_centroids
    else:
        return ('get_centroids done for ' + pdb_id)


def get_pairwise_dist(df_centroids, aa_int=None):
    """
    This function calculates pairwise distances from centroid coordinates
    that were generated from get_centroids().
    Use pdb_id to call the correct centroids.csv output file
    Returns a pandas df of pairwise distances with columns/index as AA pos
    
    df_centroids: pandas dataframe containing 4 columns of
        ['aa_num', 'x', 'y', 'z'].
    aa_int: tuple of (aa_min, aa_max) defining the aas to calculate pdists
        the default is None, which takes the min/max of df_centroids
    """

    # check for correct columns in df_centroids
    list_cols = df_centroids.columns.tolist()
    if not all(col in list_cols for col in ['aa_num', 'x', 'y', 'z']):
        raise Exception('df_centroids is missing an essential column id')
    # make sure aa_num is the correct dtype (pymol uses strings)
    df_centroids['aa_num'] = df_centroids['aa_num'].astype('int64')

    # isolate the desired amino acid interval
    # default is all resolved residues in df_centroids
    if aa_int is None:
        # remove unresolved residues (xyz = NaN) before finding aa min/max
        df_aaint = df_centroids.loc[~df_centroids.isnull().any(axis=1)].copy()
        aa_min = df_aaint['aa_num'].min()
        aa_max = df_aaint['aa_num'].max()
    else:
        aa_min = aa_int[0]
        aa_max = aa_int[1]
        df_aaint = df_centroids.loc[df_centroids['aa_num'].between(aa_min, aa_max)].copy()
        df_aaint = df_aaint.loc[~df_aaint.isnull().any(axis=1)].copy()
        if df_aaint['aa_num'].min() != aa_min:
            print('Warning! User aa_min input was ' + str(aa_min))
            print('But first resolved AA was ' + str(df_aaint['aa_num'].min()))
        if df_aaint['aa_num'].max() != aa_max:
            print('Warning! User aa_max input was ' + str(aa_max))
            print('But last resolved AA was ' + str(df_aaint['aa_num'].max()))
    # calculate all pairwise distances in euclidean 3d space
    pairwise = sp_sd.pdist(df_aaint[['x','y','z']], 'euclidean')
    # turn condensed matrix into square-form matrix
    pairwise = sp_sd.squareform(pairwise)
    # convert to pandas df with index/col as aa numbers
    df_pwdist = pd.DataFrame(pairwise, index=df_aaint['aa_num'], columns=df_aaint['aa_num'])
    return df_pwdist


def gauss(distance, std):
    arg = -(distance * distance) / (2 * std * std)
    dist = np.exp(arg)
    return dist