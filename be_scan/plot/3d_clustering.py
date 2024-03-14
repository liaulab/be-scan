# -*- coding: utf-8 -*-
"""
Created on Wed Dec  1 16:52:05 2021
@author: kevin

Script to hold the code for clustering analysis and plotting (e.g. CHUNKS, 1D clusters)

v1.0 -- 01DEC21:
    Initial script. Putting the 1D clustering stuff from David in here, and going
    to try and start re-doing the CHUNKS stuff.
v1.1 -- 03DEC21:
    Starting CHUNKS 3D analysis again. Have to try and recapitulate clusters 1/10.
v1.2 -- 10DEC21:
    Finishing optimization/parameter selection for 3D clustering for AA and gRNAs
v1.3 -- 12DEC21:
    Finalized AA and gRNA clustering parameters.
v1.4 -- 20DEC21:
    Cleaned up some sections and making plots for the clustering stuff
v2.0 -- 24FEB22:
    Adjusting size of Fig 2a (1D clustering plot)
v2.1 -- 08APR22:
    Making a new PyMOL panel for Figure 2, don't actually know why I updated this version
v2.2 -- 26APR22:
    Making figure panels for Supplemental Figure 2. Also removed a lot of scratch code.
    See previous versions for missing scratch code.
"""
#%% import packages

import os, re, sys, time, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.stats as stats
import scipy.spatial.distance as dist
import scipy.interpolate as interp
import scipy.optimize as sp_opt
import scipy.cluster as sp_cl


# import statsmodels.api as sm
import statsmodels.stats.multitest as smm
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from Bio.Data import IUPACData

sns.set_context('talk')
#%% internal functions for analysis (taken from KN10164_analysis_v2_211127.py)

# fx to avg PROVEAN scores for gRNAs cutting between codons -- v2, updated 25JAN21
def avg_prov(gene, cut_pos, df_scores, aa_col='aa_pos', score_col='provean'):
    df_scores = df_scores.loc[df_scores['Gene'] == gene].copy()
    aa_max = df_scores.loc[df_scores['Gene'] == gene][aa_col].max()
    if cut_pos > (aa_max + 1):
        raise Exception('cut_pos is greater than aa_max + 1')
    elif cut_pos == (aa_max + 1): # if it hits stop codon, then assign last AA value
        cut_pos = aa_max
    # if cut_pos hits single AA, assign score; else, avg AA on either side of cut_pos
    if cut_pos % 1 == 0:
        score = df_scores.loc[df_scores[aa_col] == cut_pos, score_col].values[0]
    else:
        aa1, aa2 = cut_pos - 0.5, cut_pos + 0.5
        if cut_pos == 0.5:
            score = df_scores.loc[df_scores[aa_col] == 1, score_col].values[0]
        elif cut_pos == aa_max + 0.5:
            score = df_scores.loc[df_scores[aa_col] == aa_max, score_col].values[0]
        else:
            score = df_scores.loc[(df_scores[aa_col] == aa1) | (df_scores[aa_col] == aa2)][score_col].mean()
    return score

# revised loess fx v3; allow choice between statsmodels interp or scipy interp1d
def loess_v3(x_obs, y_obs, span, x_out=None, interp_how='quadratic', it=0, loess_kws=dict(missing='raise',return_sorted=False), interp_kws=dict(fill_value='extrapolate')):
    x_obs = x_obs.astype(float)
    if x_out is None:
        x_out = x_obs
    df_loess = pd.DataFrame()
    df_loess['xvals'] = x_out.astype(float)
    df_loess = df_loess.sort_values(by='xvals').reset_index(drop=True)
    if interp_how == 'statsmodels':
        df_loess['y_loess'] = lowess(endog=y_obs, exog=x_obs, xvals=df_loess['xvals'], frac=span, it=it, **loess_kws)
    else:
        df_interp = pd.DataFrame()
        df_interp['x_obs'] = x_obs
        df_interp['y_obs_loess'] = lowess(endog=y_obs, exog=x_obs, xvals=x_obs, frac=span, it=it, **loess_kws)
        df_interp = df_interp.groupby('x_obs',as_index=False).agg({'y_obs_loess':'mean'})
        fx_interp = interp.interp1d(x=df_interp['x_obs'], y=df_interp['y_obs_loess'], kind=interp_how, **interp_kws)
        df_loess['y_loess'] = fx_interp(df_loess['xvals'])
    df_loess['type'] = np.where(df_loess['xvals'].isin(x_obs), 'loess', 'interp')
    return df_loess

# fx to generate gRNA annotations for plotting
def label_grnas(cut_pos, gene, df_ref):
    if cut_pos % 1 != 0.5:
        resi = df_ref.loc[(df_ref['Gene'] == gene) & (df_ref['aa_pos'] == cut_pos), 'aa_ref'].values[0]
        label = 'g' + resi + str(int(cut_pos))
    else:
        aa1, aa2 = int(cut_pos - 0.5), int(cut_pos + 0.5)
        # edge cases where cut_site_AA == 0.5 or max AA.5 (b/w last codon and stop)
        if aa1 == 0:
            resi = df_ref.loc[(df_ref['Gene'] == gene) & (df_ref['aa_pos'] == aa2), 'aa_ref'].values[0]
            label = 'g' + resi + str(aa2)
        elif aa2 > df_ref.loc[df_ref['Gene'] == gene, 'aa_pos'].max():
            resi = df_ref.loc[(df_ref['Gene'] == gene) & (df_ref['aa_pos'] == aa1), 'aa_ref'].values[0]
            label = 'g' + resi + str(aa1)
        else:
            resi1 = df_ref.loc[(df_ref['Gene'] == gene) & (df_ref['aa_pos'] == aa1), 'aa_ref'].values[0]
            resi2 = df_ref.loc[(df_ref['Gene'] == gene) & (df_ref['aa_pos'] == aa2), 'aa_ref'].values[0]
            label = 'g' + resi1 + str(aa1) + '/' + resi2 + str(aa2)
    return label

# fx to color the domain intervals on CDS scatterplots
def color_domains(ax, domains, palette, spines=True, dom_kws=dict(alpha=0.3), nodom_kws=dict(fc='#CCCCCC', alpha=0.3)):
    for (xmin,xmax),color in zip(list(domains.values())[:-1],list(palette.values())[:-1]):
        ax.axvspan(xmin=xmin-0.5, xmax=xmax+0.5, fc=color, **dom_kws)
    for xmin,xmax in domains['Undefined']:
        ax.axvspan(xmin=xmin-0.5, xmax=xmax+0.5, **nodom_kws)
    if spines:
        for edge,spine in ax.spines.items():
            spine.set_visible(True)
    return ax

# fx to annotate domains in the input data
def annotate_domains(df_input, aa_col, dict_domains):
    df = pd.DataFrame()
    df['aa_pos'] = df_input[aa_col]
    df['Domain'] = np.nan
    for domain,xvals in dict_domains.items():
        if domain != 'Undefined':
            xmin, xmax = xvals[0], xvals[1]
            df['Domain'] = np.where(df['aa_pos'].between(xmin, xmax), domain, df['Domain'])
        elif domain == 'Undefined':
            for xmin,xmax in xvals:
                df['Domain']  = np.where(df['aa_pos'].between(xmin, xmax), 'Undefined', df['Domain'])
    if df['Domain'].isna().any():
        print('Warning! NaNs found in output "Domain" column')
    return df['Domain']

### fx to adjust widths of seaborn boxplot boxes (from github)
### https://github.com/mwaskom/seaborn/issues/1076
### 08DEC21: v2 -- adapted to allow calling axes instead of fig

from matplotlib.patches import PathPatch
def adjust_box_widths(g, fac, is_ax=False):
    """
    Adjust the widths of a seaborn-generated boxplot.
    g: matplotlib figure or axes object
        if passing an axes object, be sure to call is_ax=True
    fac: scaling factor (fraction)
    is_ax: bool, default False
        whether the passed object (g) is a figure or axes object
    """
    if is_ax: # if g is an axes object
        ax = g
        # iterating through axes artists:
        for c in ax.get_children():
            # searching for PathPatches
            if isinstance(c, PathPatch):
                # getting current width of box:
                p = c.get_path()
                verts = p.vertices
                verts_sub = verts[:-1]
                xmin = np.min(verts_sub[:, 0])
                xmax = np.max(verts_sub[:, 0])
                xmid = 0.5*(xmin+xmax)
                xhalf = 0.5*(xmax - xmin)
                # setting new width of box
                xmin_new = xmid-fac*xhalf
                xmax_new = xmid+fac*xhalf
                verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new
                # setting new width of median line
                for l in ax.lines:
                    if np.all(l.get_xdata() == [xmin, xmax]):
                        l.set_xdata([xmin_new, xmax_new])
    else: # if g is a figure object
        # iterating through Axes instances
        for ax in g.axes:
            # iterating through axes artists:
            for c in ax.get_children():
                # searching for PathPatches
                if isinstance(c, PathPatch):
                    # getting current width of box:
                    p = c.get_path()
                    verts = p.vertices
                    verts_sub = verts[:-1]
                    xmin = np.min(verts_sub[:, 0])
                    xmax = np.max(verts_sub[:, 0])
                    xmid = 0.5*(xmin+xmax)
                    xhalf = 0.5*(xmax - xmin)
                    # setting new width of box
                    xmin_new = xmid-fac*xhalf
                    xmax_new = xmid+fac*xhalf
                    verts_sub[verts_sub[:, 0] == xmin, 0] = xmin_new
                    verts_sub[verts_sub[:, 0] == xmax, 0] = xmax_new
                    # setting new width of median line
                    for l in ax.lines:
                        if np.all(l.get_xdata() == [xmin, xmax]):
                            l.set_xdata([xmin_new, xmax_new])

#%% import pre-processed data (processed in KN10164_analysis_v2_211127.py)

# subdirectory to hold analysis/output files (05OCT21)
outpath = Path.cwd() /'KN10164_analysis_v1_211004'
# defining lists of common column headers for downstream ease
list_refcols = ['sgRNA_ID', 'sgRNA_seq', 'Gene', 'cut_site_NT', 'cut_site_AA', 'Domain', 'spec_score', 'doench_score']
list_timepoints = ['7d','14d','21d','28d','42d','56d']
list_lfccols = [y + x for y in ['DAC_','DMSO_'] for x in list_timepoints]
list_enrichcols = ['enrich_7d', 'enrich_14d', 'enrich_21d', 'enrich_28d', 'enrich_42d', 'enrich_56d']
list_condcols = ['fit_7d', 'fit_14d', 'fit_21d', 'fit_28d', 'fit_42d', 'fit_56d', 'func_7d', 'func_14d', 'func_21d', 'func_28d', 'func_42d', 'func_56d']

### importing processed CRISPR screen data for DNMT1/UHRF1 ###

### 29NOV21: import v2 data (fixed enriched gRNA calling)
# master screen data (all scores, incl. pre-norm and post-norm)
in_lfc = pd.read_csv('./KN10164_analysis_v1_211004/KN10164_processed_data/' + 'KN10164_screen_master_v2_211129.csv')
# selected processed screen data (ctrl-norm, dmso-norm (func) scores)
df_lfc = pd.read_csv('./KN10164_analysis_v1_211004/KN10164_processed_data/' + 'KN10164_screen_norm_v2_211129.csv')
# annotated and processed data, no ctrl gRNAs
df_dnmt = pd.read_csv('./KN10164_analysis_v1_211004/KN10164_processed_data/' + 'KN10164_screen_norm_dnmt1_v2_211129.csv').reset_index(drop=True)
df_uhrf = pd.read_csv('./KN10164_analysis_v1_211004/KN10164_processed_data/' + 'KN10164_screen_norm_uhrf1_v2_211129.csv').reset_index(drop=True)
# negative control gRNA statistics (e.g. mean, mean + 2sd)
df_ctrl = pd.read_csv('./KN10164_analysis_v1_211004/KN10164_processed_data/' + 'KN10164_screen_ctrl_stats_v1_211005.csv')

### importing PROVEAN/ConSurf/gnomAD data for DNMT1/UHRF1 ###
# aggregated, processed, and consolidated provean/consurf/gnomad data
df_cons = pd.read_csv('./KN10164_analysis_v1_211004/KN10164_processed_data/' + 'KN10164_agg_conservation_scores_v1_211005.csv')

#%% setting universal plot parameters

### setting inches to mm conversion ratio (mpl figsize in inches, AI in mm)
mm = 1/25.4
### setting seaborn context to 'paper' and making matplotlib rcParams adjustments
plot_rc = {'font.size': 6.0, 'axes.titlesize': 7.0, 'axes.labelsize': 6.0, 'xtick.labelsize': 5.0, ### FONT SIZES
           'ytick.labelsize': 5.0,'legend.fontsize': 6.0, 'legend.title_fontsize': 7.0, ### FONT SIZES
           'axes.linewidth': 0.5, 'grid.linewidth': 0.4, 'lines.linewidth': 1.0, 'lines.markersize': 2.5, ### LINE/MARKER SIZES
           'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2.4, 'ytick.major.size': 2.4} ### LINE/MARKER SIZES
sns.set_context('paper', rc=plot_rc)
### setting other matplotlib rcParams (e.g. text/label padding)
mpl.rcParams['axes.labelpad'] = 2.0
mpl.rcParams['xtick.major.pad'] = 0.85
mpl.rcParams['ytick.major.pad'] = 0.85

### in case you need to reset to 'talk' context (see mpl style sheet)
# sns.set_context('talk')
# mpl.style.use('KCN_figures_talk_v1_211201')

### 29NOV21: removing bbox_inches='tight' because this affects figure sizing
plot_kws_savepng = dict(format='png', dpi=300, transparent=False)
plot_kws_savepdf = dict(format='pdf', transparent=True)

#%% 03DEC21: internal fx for CHUNKS 3D clustering

# get_pairwise_dist fx from KCN_CHUNKS_clusteringfunctions_v3_201104.py
def get_pairwise_dist(df_centroids, aa_int=None):
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
    pairwise = dist.pdist(df_aaint[['x','y','z']], 'euclidean')
    # turn condensed matrix into square-form matrix
    pairwise = dist.squareform(pairwise)
    # convert to pandas df with index/col as aa numbers
    df_pwdist = pd.DataFrame(pairwise, index=df_aaint['aa_num'], columns=df_aaint['aa_num'])
    return df_pwdist

# hill fx for scaling LFCs from KCN_CHUNKS_clusteringfunctions_v3_201104.py
def hill(lfc, m, theta):
    num = lfc**m
    denom = (lfc**m) + (theta**m)
    val = num/denom
    return val

# gaussian fx for scaling LFCs from KCN_CHUNKS_clusteringfunctions_v3_201104.py
def gauss(distance, std):
    arg = -(distance * distance) / (2 * std * std)
    dist = np.exp(arg)
    return dist

# rounding sgRNAs with cut site between codons from KCN_CHUNKS_clusteringfunctions_v3_201104.py
def chunks_round(df_col, how):
    if how == 'default':
        new_col = df_col.round()
    elif how == 'up':
        new_col = df_col.apply(np.ceil)
    elif how == 'down':
        new_col = df_col.apply(np.floor)
    else:
        raise Exception('Invalid value for how. Must be default, up or down')
    return new_col

# CHUNKS scaling fx from KCN_CHUNKS_clusteringfunctions_v3_201104.py
def chunks_scale(in_lfc, df_pwdist, cond, gene, sg_round='default', std=16, m=2, theta=3):
    # import the log2 fold change data if path or str
    if isinstance(in_lfc, pd.DataFrame):
        df_import = in_lfc.copy()
    else:
        df_import = pd.read_csv(in_lfc)
    # check for essential columns
    list_reqcols = ['Gene', 'cut_site_AA', 'sgRNA_ID', cond]
    if not all(col in df_import.columns.tolist() for col in list_reqcols):
        list_miss = [col for col in list_reqcols if col not in df_import.columns.tolist()]
        raise Exception('in_lfc missing essential column(s): ' + str(list_miss))
    # isolate enrichment values for the gene of interest
    df_lfc = df_import.loc[df_import['Gene'] == gene].copy()
    # decide how to deal with sgRNAs that cut between codons (tweeners)
    # default is to round up, alternative is to round down
    df_lfc['aa_num'] = chunks_round(df_col=df_lfc['cut_site_AA'], how=sg_round)
    # remove all sgRNAs that target AAs without pairwise distances
    df_lfc = df_lfc.loc[df_lfc['aa_num'].isin(df_pwdist.index)].copy()
    # sort by AA position, then re-index
    df_lfc = df_lfc.sort_values(by=['aa_num', 'sgRNA_ID']).reset_index()
    ids = df_lfc['sgRNA_ID']

    # generate a pairwise matrix of summed LFC values with matrix operation
    df_sum = df_lfc[cond].values[:, None] + df_lfc[cond].values[None, :]
    df_sum = pd.DataFrame(index=ids, columns=ids, data=df_sum)
    # generate a matrix to remember the sign
    df_sign = df_sum.apply(lambda x: np.where(x > 0, 1, -1))
    # get abs value, scale w/ hill function, then adjust sign
    df_sum = df_sum.abs()
    df_sum = df_sum.apply(lambda x: hill(x, m, theta))
    df_hill = df_sum * df_sign
    # get AAs using sgRNA ids. not the best way but i'm overly anxious about the order
    aas = df_lfc.loc[df_lfc['sgRNA_ID'] == ids]['aa_num']
    # retrieve pairwise distances from df_pwdist for the relevant sgRNAs
    df_gauss = df_pwdist.loc[aas, aas].copy()
    # reset index/column to match sgRNA ID rather than AA number
    df_gauss.set_index(keys=ids, drop=True, inplace=True)
    df_gauss.set_axis(labels=ids, axis=1, inplace=True)
    # scale the pairwise distances with the gaussian function
    df_gauss = df_gauss.apply(lambda x: gauss(x, std))
    # multiply gaussian component w/ hill component to get scaled LFCs
    df_scaled = df_gauss * df_hill
    # change data type to float
    df_scaled = df_scaled.astype(float)

    return df_scaled

#%% 03DEC21: import CHUNKS 3D clustering data/define variables

### calculate pairwise distances from KN10152_4wxx_centroids_v2.csv (1st time only, then import)
# df_cent = pd.read_csv('./KN10164_clustering_analysis_v1_211201/KN10152_4wxx_centroids_v2.csv')
# df_pwdist = get_pairwise_dist(df_cent)
# df_pwdist.to_csv('./KN10164_clustering_analysis_v1_211201/KN10164_4wxx_pairwise_dist_v1_211203.csv')
df_pwdist = pd.read_csv('./KN10164_clustering_analysis_v1_211201/KN10164_4wxx_pairwise_dist_v1_211203.csv', index_col=0)
df_pwdist.columns = df_pwdist.columns.astype(int) # make columns int64 dtype

### get list of gRNAs tested in BL106
list_bl106 = ['DNMT1_' + str(x) for x in [197,205,220,229,247,251,287,310,472,477,514,515,529,581,594,600,620,624,627,644,650,679,709,722,746,752,768,808]]

### get distance scaling factors by scaling pw dist with gaussian fx
# using std=16 angstrom for gRNAs, using std=10 for AA LOESS
df_gauss = df_pwdist.apply(lambda x: gauss(x, std=16)).copy()
# df_gauss2 = df_pwdist.apply(lambda x: gauss(x, std=10)).copy()

### import the pre-calculated CHUNKS clustering results for gRNAs and AAs
df_clus = pd.read_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_grna_clusters_v1_211213.csv')
df_clstats = pd.read_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_grna_clusters_stats_v1_211213.csv')

df_clusls = pd.read_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_aa_loess_clusters_v1_211213.csv')
df_clsagg = pd.read_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_aa_loess_clusters_stats_v1_211213.csv')

### get list of gRNA target AAs in new CHUNKS clusters 1 and 2
# list_cl1 = [511,528,531,545,546,557,560,566,567,576,585,587,588,593,594,595,649,650,652,693,696,697,702,1425,1479,1489,1503,1504,1505,1532]
# list_cl2 = [1077,1081,1082,1143,1147,1152,1157,1158,1159,1160,1266,1299,1303,1326,1327,1558,1587]

#%% 26APR22: CHUNKS analysis for Supp. Fig. 2 -- PWES Hi-C heatmap


###### ANALYSIS 1: Perform CHUNKS analysis to get PWES score matrix for all resolved sgRNAs ######
    #### SUMMARY WORKFLOW:
    #### get pairwise distances for 4wxx centroids from KN10152 (df_pwdist)
    #### calculate gaussian pairwise distance matrix with std=16 (df_gauss)
    #### get pairwise sums of resistance scores for ALL DNMT1 sgRNAs (incl. non-resolved sgRNAs)
    #### z-score the pairwise sums using ALL DNMT1 sgRNAs
    #### scale the z-scores with tanh fx and then scale with gaussian fx (df_pws1)


    ### import pre-calculated pairwise distances (calculated from KN10152_4wxx_centroids_v2.csv)
df_pwdist = pd.read_csv('./KN10164_clustering_analysis_v1_211201/KN10164_4wxx_pairwise_dist_v1_211203.csv', index_col=0)
df_pwdist.columns = df_pwdist.columns.astype(int) # make columns int64 dtype

    ### scale pairwise distances with gaussian fx (std=16 angstroms) ###
df_gauss = df_pwdist.apply(lambda x: gauss(x, std=16)).copy()

    ### calculate pairwise summed scores for sgRNAs ###
# first calculate pairwise sums of resistance scores for ALL DNMT1 sgRNAs (incl. non-resolved sgRNAs)
df_scores = df_dnmt.copy()
df_scores['aa_pos'] = df_scores['cut_site_AA'].round() # round grna cut sites
list_aas = df_scores[df_scores['aa_pos'].isin(df_gauss.index)]['aa_pos'] # get all resolved gRNAs (n=646)
# make pairwise sums matrix for resolved sgRNAs
df_pws = df_scores[df_scores['aa_pos'].isin(list_aas)].copy() # get resolved gRNAs (n=646)
df_pws = pd.DataFrame(index=df_pws['sgRNA_ID'], columns=df_pws['sgRNA_ID'],
                      data=df_pws['func_56d'].values[:, None] + df_pws['func_56d'].values[None, :]) # pairwise sums matrix
# calculate pairwise sums for ALL DNMT1 sgRNAs (n=830) to find mean and stdev for z-scoring and tanh-scaling
temp = pd.DataFrame(index=df_scores['sgRNA_ID'], columns=df_scores['sgRNA_ID'], data=df_scores['func_56d'].values[:, None] + df_scores['func_56d'].values[None, :])
temp = pd.DataFrame(index=temp.index, columns=temp.columns, data=np.where(np.triu(np.ones(temp.shape), k=1).astype(bool), temp, np.nan))
temp2 = pd.Series([y for x in temp.columns for y in temp[x]], name='sum_lfc').dropna() # pw summed gRNA func scores
# z-score resolved sgRNA scores (df_pws) w/ pairwise sum mean and std of all sgRNAs, then tanh-scale the z-scores
df_pws1 = np.tanh((df_pws - temp2.mean()) / temp2.std()) # all gRNAs (n=830) pw_sum mean == -0.845; std == 4.406
df_pws1.index, df_pws1.columns = list_aas, list_aas # replace index/columns
df_pws1 = df_pws1 * df_gauss.loc[list_aas, list_aas].copy()


#%% 26APR22: CHUNKS analysis for Supp. Fig. 2 -- WAP and cliuster significance

###### ANALYSIS 2: Simulate null distribution of WAPs by shuffling AA positions summed PWES ######
    #### SUMMARY WORKFLOW:
    #### shuffle sgRNA aa positions and re-calculate PWES
    #### calculate the WAP using the absolute value of the sgRNAs per cluster

    ### generate shuffled aa_positions for resolved sgRNAs (n=10000)
df_sim = df_scores[df_scores['aa_pos'].isin(df_gauss.index)][['sgRNA_ID','func_56d','aa_pos']].copy() # get all resolved sgRNAs (n=646)
df_sim = df_sim.rename(columns={'aa_pos':'aa_obs'}).reset_index(drop=True) 
temp_dict = {}
for i in np.arange(0,10000):
    temp_dict['iter' + str(i)] = df_sim['aa_obs'].sample(frac=1, ignore_index=True).rename('iter' + str(i))
df_sim = pd.concat([df_sim] + [x for x in temp_dict.values()], axis=1)
# df_sim.to_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_shuffled_aa_pos_v1_220428.csv', index=False)


    ### set various "universal" variables
### calculate the pairwise sums matrix, z-score, tanh-transform
# make pairwise sums matrix for resolved sgRNAs
temp_pws = pd.DataFrame(index=df_sim['sgRNA_ID'], columns=df_sim['sgRNA_ID'],
                        data=df_sim['func_56d'].values[:, None] + df_sim['func_56d'].values[None, :])
# calculate pairwise sums for ALL DNMT1 sgRNAs (n=830) to find mean and stdev for z-scoring and tanh-scaling
temp = pd.DataFrame(index=df_scores['sgRNA_ID'], columns=df_scores['sgRNA_ID'], data=df_scores['func_56d'].values[:, None] + df_scores['func_56d'].values[None, :])
temp = pd.DataFrame(index=temp.index, columns=temp.columns, data=np.where(np.triu(np.ones(temp.shape), k=1).astype(bool), temp, np.nan))
temp2 = pd.Series([y for x in temp.columns for y in temp[x]], name='sum_lfc').dropna() # pw summed gRNA func scores
# z-score resolved sgRNA scores (df_pws) w/ pairwise sum mean and std of all sgRNAs, then tanh-scale the z-scores
temp_pws = np.tanh((temp_pws - temp2.mean()) / temp2.std()) # all gRNAs (n=830) pw_sum mean == -0.845; std == 4.406
### create a 646x646 mask to set self vs. self diagonal to 0
temp_mask = np.tril(np.ones(temp_pws.shape)).astype(bool)
### make dict of (cluster, list of sgRNAs)
dict_clus = {k:df_clus[df_clus['cl_new_rank'] == k]['sgRNA_ID'].tolist() for k in np.arange(1,20)}


    ### calculate the observed WAP for each cluster
# create new gaussian matrix, scale the pw sum matrix, get abs value
temp_gauss = df_gauss.loc[df_sim['aa_obs'], df_sim['aa_obs']].copy()
# reset index and change to sgRNA_ID index/columns
temp_gauss.index, temp_gauss.columns = np.arange(0,646), np.arange(0,646)
temp_gauss = temp_gauss.rename(index=df_sim['sgRNA_ID'], columns=df_sim['sgRNA_ID'])
# scale matrix, get absolute value, mask lower triangle
temp_pws2 = (temp_pws * temp_gauss).abs()
temp_pws2 = temp_pws2.mask(temp_mask, 0)
# calculate the observed WAP for each cluster -- use only within cluster PWES interactions
df_wap = pd.DataFrame(columns=['cluster', 'obs_wap'])
df_wap['cluster'] = np.arange(1,20)
df_wap['obs_wap'] = pd.Series([temp_pws2.loc[dict_clus[x],dict_clus[x]].sum().sum() for x in dict_clus])


    ### calculate the simulated WAP distribution for each cluster (n=10000)
temp_dict = {}
# for each iteration, follow same workflow as for the observed WAP
for i in np.arange(0,10000):
    temp_gauss = df_gauss.loc[df_sim['iter' + str(i)], df_sim['iter' + str(i)]].copy()
    temp_gauss.index, temp_gauss.columns = np.arange(0,646), np.arange(0,646)
    temp_gauss = temp_gauss.rename(index=df_sim['sgRNA_ID'], columns=df_sim['sgRNA_ID'])
    temp_pws2 = (temp_pws * temp_gauss).abs()
    temp_pws2 = temp_pws2.mask(temp_mask, 0)
    temp_dict['iter' + str(i)] = pd.Series([temp_pws2.loc[dict_clus[x],dict_clus[x]].sum().sum() for x in dict_clus], name='iter' + str(i))
df_wap = pd.concat([df_wap] + [x for x in temp_dict.values()], axis=1)
# df_wap.to_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_shuffled_aa_pos_null_waps_v1_220428.csv', index=False)

    ### calculate stats/p-vals for observed WAP vs. null distribution
df_wapst = df_wap[['cluster','obs_wap']].copy()
df_wapst['mean'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].mean(axis=1)
df_wapst['std'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].std(axis=1)
df_wapst['95ci_min'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].quantile(q=0.025, axis=1)
df_wapst['95ci_max'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].quantile(q=0.975, axis=1)
df_wapst['obs_gt'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].gt(df_wap['obs_wap'], axis=0).sum(axis=1)
df_wapst['1t_pval'] = df_wapst['obs_gt'] / 10000
# df_wapst.to_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_shuffled_aa_pos_wapstats_v1_220428.csv', index=False)


#%% 26APR22: Making plots for Supplemental Figure 2 (formatted for AI)

### Supplemental Fig. 2a -- Hi-C heatmap of the sgRNA PWES matrix
plot_kws_heatmap = dict(vmin=-1, vmax=1, center=0, cmap='RdBu_r', square=True, rasterized=True, clip_on=False,
                        xticklabels=False, yticklabels=False,
                        cbar_kws=dict(ticks=[-1,0,1], label='PWES', orientation='vertical', fraction=0.1, shrink=0.6, aspect=8, clip_on=False))
plot_mask = np.zeros_like(df_pws1)
plot_mask[np.tril_indices_from(df_pws1)] = True
# draw domain cutoffs with vertical lines on x axis
temp = df_cons[df_cons['Gene'] == 'DNMT1'].copy().reset_index(drop=True)
temp2 = pd.DataFrame()
temp2['aa_pos'] = df_pws1.index.copy()
temp2['Domain'] = temp2['aa_pos'].apply(lambda x: temp[temp['aa_pos'] == x]['Domain'].values[0])
temp3 = {}
for plot_dom in temp2['Domain'].unique():
    if plot_dom != 'Undefined':
        temp3[plot_dom] = (temp2[temp2['Domain'] == plot_dom]['aa_pos'].idxmin(), temp2[temp2['Domain'] == plot_dom]['aa_pos'].idxmax())
# Supplemental Fig. 2a plotting
fig, ax = plt.subplots(figsize=(65*mm,50*mm))
sns.heatmap(ax=ax, data=df_pws1, **plot_kws_heatmap)#, mask=plot_mask)
ax.xaxis.set_visible(False), ax.yaxis.set_visible(False)
for plot_spine in ax.spines.values():
    plot_spine.set(color='k', visible=True)
ax._colorbars[0].spines['outline'].set(color='k', linewidth=0.5, visible=True) # set colorbar spines
ax.axline((0,0), slope=1, c='k', lw=0.5) # draw diagonal line
for plot_dom, (plot_min, plot_max) in temp3.items(): # draw domain cutoffs
    if plot_min != 0:
        ax.axvline(x=plot_min, c=pal_dnmt[plot_dom], lw=0.5, clip_on=False)
    if plot_max != 646:
        ax.axvline(x=plot_max + 1, c=pal_dnmt[plot_dom], lw=0.5, clip_on=False)
# fig.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Supp_Fig2_3d_clustering_hi-c_heatmap_v2_220426.pdf', **plot_kws_savepdf)


### Supplemental Fig. 2b -- obs vs. simulated WAPs histogram/CDFs
temp_plot = df_wap.transpose()
temp_plot.columns = temp_plot.loc['cluster']
temp_plot = temp_plot.drop(index='cluster')
# plotting
plot_kws_hist = dict(stat='frequency', bins=100, kde=False, fill=True, color='#999999', edgecolor='none', alpha=0.8, clip_on=False)
fig, ((ax1,ax2), (ax3,ax4)) = plt.subplots(ncols=2, nrows=2, figsize=(100*mm,100*mm))
sns.histplot(ax=ax1, data=temp_plot.loc[['iter' + str(i) for i in np.arange(0,10000)]], x=1, binrange=(0,200), **plot_kws_hist)
sns.ecdfplot(ax=ax2, data=temp_plot.loc[['iter' + str(i) for i in np.arange(0,10000)]], x=1, clip_on=False, color='#525252')
sns.histplot(ax=ax3, data=temp_plot.loc[['iter' + str(i) for i in np.arange(0,10000)]], x=2, binrange=(0,70), **plot_kws_hist)
sns.ecdfplot(ax=ax4, data=temp_plot.loc[['iter' + str(i) for i in np.arange(0,10000)]], x=2, clip_on=False, color='#525252')
for ax in ax1,ax2,ax3,ax4:
    ax.grid(visible=False)
    ax.xaxis.set_clip_on(False), ax.yaxis.set_clip_on(False)
    for plot_spine in ax.spines.values():
        plot_spine.set(color='k', visible=True)
ax1.axvline(x=temp_plot.at['obs_wap', 1], c='#e41a1c', lw=1)
ax1.set(title='Cluster 1 sgRNAs', xlabel='Summed PWES', ylabel='Frequency', xlim=(0,200), ylim=(0,600))
ax2.set(xlabel='Summed PWES', ylabel='Cumulative probability', xlim=(20,100), ylim=(-0.02,1.02))
ax3.axvline(x=temp_plot.at['obs_wap', 2], c='#e41a1c', lw=1)
ax3.set(title='Cluster 2 sgRNAs', xlabel='Summed PWES', ylabel='Frequency', xlim=(0,70), ylim=(0,1400))
ax4.set(xlabel='Summed PWES', ylabel='Cumulative probability', xlim=(0,40), ylim=(-0.02,1.02))
# fig.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Supp_Fig2_3d_clus1_clus2_wap_hist_ecdf_v1_220428.pdf', **plot_kws_savepdf)


#%% 26APR22: CHUNKS analysis scratch


###### NOTE: THIS IS INCORRECT AND WE SHOULD SHUFFLE THE AA POSITIONS NOT THE SGRNAs
###### ANALYSIS 2: Simulate CHUNKS and get null distribution of summed PWES ######
    #### SUMMARY WORKFLOW:
    #### shuffle sgRNA resistance scores (RESOLVED sgRNAs ONLY) and calculate pairwise sums matrix
    #### z-score pairwise sums using mean/std of all sgRNAs
    #### tanh-transform z-scores and then scale with gaussian
    #### calculate the WAP using the absolute value of the sgRNAs per cluster


    ### generate shuffled resistance scores for resolved sgRNAs (n=1000)
df_sim = df_scores[df_scores['aa_pos'].isin(df_gauss.index)][['sgRNA_ID','aa_pos','func_56d']].copy()
df_sim = df_sim.rename(columns={'func_56d':'obs'}).reset_index(drop=True) # get all resolved sgRNAs (n=646)
temp_dict = {}
for i in np.arange(0,10000):
    temp_dict['iter' + str(i)] = df_sim['obs'].sample(frac=1, ignore_index=True).rename('iter' + str(i))
df_sim = pd.concat([df_sim] + [x for x in temp_dict.values()], axis=1)
# df_sim.to_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_shuffled_scores_null_waps_v1_220427.csv', index=False)

    ### set various "universal" variables
### calculate pairwise sums for ALL DNMT1 sgRNAs (n=830) to find mean and stdev for z-scoring and tanh-scaling
temp = pd.DataFrame(index=df_scores['sgRNA_ID'], columns=df_scores['sgRNA_ID'], data=df_scores['func_56d'].values[:, None] + df_scores['func_56d'].values[None, :])
temp = pd.DataFrame(index=temp.index, columns=temp.columns, data=np.where(np.triu(np.ones(temp.shape), k=1).astype(bool), temp, np.nan))
temp2 = pd.Series([y for x in temp.columns for y in temp[x]], name='sum_lfc').dropna() # pw summed gRNA func scores
### create 646x646 dataframe of the pairwise gaussian scaling factors for the resolved sgRNAs
# rename the index/columns to the sgRNA ID for easier selection later on
df_gres = df_gauss.loc[list_aas, list_aas].copy()
df_gres.index, df_gres.columns = np.arange(0,646), np.arange(0,646) # need to reset index first bc of duplicate AAs
df_gres = df_gres.rename(index=df_sim['sgRNA_ID'], columns=df_sim['sgRNA_ID'])
### create a 646x646 mask to set self vs. self diagonal to 0
temp_mask = np.zeros_like(df_gres)
temp_mask[np.diag_indices_from(temp_mask)] = 1
temp_mask = np.where(temp_mask == 1, True, False)
### make dict of (cluster, list of sgRNAs)
dict_clus = {k:df_clus[df_clus['cl_new_rank'] == k]['sgRNA_ID'].tolist() for k in np.arange(1,20)}


    ### calculate the observed WAP for each cluster
# create pw sum matrix, z-score and tanh scale, scale by gaussian, convert to absolute value
temp_pws = pd.DataFrame(index=df_sim['sgRNA_ID'], columns=df_sim['sgRNA_ID'],
                        data=df_sim['obs'].values[:, None] + df_sim['obs'].values[None, :])
temp_pws = np.tanh((temp_pws - temp2.mean()) / temp2.std())
temp_pws = (temp_pws * df_gres).copy().abs()
temp_pws = temp_pws.mask(temp_mask, 0) # set self vs. self diagonal to 0
# calculate the observed WAP for each cluster
df_wap = pd.DataFrame(columns=['cluster', 'obs_wap'])
df_wap['cluster'] = np.arange(1,20)
df_wap['obs_wap'] = pd.Series([temp_pws[dict_clus[x]].sum().sum() for x in dict_clus])

    ### calculate the simulated WAP distribution for each cluster (n=10000)
temp_dict = {}
# for each iteration, follow same workflow as for the observed WAP
for i in np.arange(0,10000):
    temp_pws = pd.DataFrame(index=df_sim['sgRNA_ID'], columns=df_sim['sgRNA_ID'],
                            data=df_sim['iter' + str(i)].values[:, None] + df_sim['iter' + str(i)].values[None, :])
    temp_pws = np.tanh((temp_pws - temp2.mean()) / temp2.std())
    temp_pws = (temp_pws * df_gres).copy().abs()
    temp_pws = temp_pws.mask(temp_mask, 0) # set self vs. self diagonal to 0
    temp_dict['iter' + str(i)] = pd.Series([temp_pws[dict_clus[x]].sum().sum() for x in dict_clus], name='iter' + str(i))
df_wap = pd.concat([df_wap] + [x for x in temp_dict.values()], axis=1)
# df_wap.to_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_shuffled_scores_null_waps_v1_220427.csv', index=False)

    ### calculate stats/p-vals for observed WAP vs. null distribution
df_wapst = df_wap[['cluster','obs_wap']].copy()
df_wapst['mean'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].mean(axis=1)
df_wapst['std'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].std(axis=1)
df_wapst['95ci_min'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].quantile(q=0.025, axis=1)
df_wapst['95ci_max'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].quantile(q=0.975, axis=1)
df_wapst['obs_gt'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].gt(df_wap['obs_wap'], axis=0).sum(axis=1)
df_wapst['1t_pval'] = df_wapst['obs_gt'] / 10000


#%% 15DEC21: Figure 2 plots (1D/3D clustering)

plot_kws_scatter = dict(color='#525252', alpha=0.3, s=9, ec='none', clip_on=False, legend=None)
plot_kws_line = dict(color='#e41a1c', alpha=1, legend=None)

### Fig 2a schematic -- DNMT1 56d function score vs. CDS scatter plot (all gray, no LOESS overlay)
fig, ax = plt.subplots(figsize=(32*mm, 24*mm)) # size is (width, height)
sns.scatterplot(ax=ax, data=df_dnmt, x='cut_site_AA', y='func_56d', zorder=1.5, **plot_kws_scatter)
sns.lineplot(ax=ax, data=dnmtls, x='aa_pos', y='func_56d_ls', zorder=2, **plot_kws_line)
ax.set(xlabel=None, xlim=(1,1616), xticks=[], ylabel=None, ylim=(-12,12), yticks=[])
ax.grid(b=None)
color_domains(ax=ax, domains=doms_dnmt, palette=pal_dnmt, spines=False)
fig.set_size_inches(32*mm, 24*mm) # sometimes need to resize due to monitor issues
# plt.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Fig2_schematic_1d_clus_func_score_scatter_v2_211215.pdf', **plot_kws_savepdf)

### Fig 2a schematic -- shuffled gRNA LOESS overlays for N=20 iterations
fig, ax = plt.subplots(figsize=(32*mm, 24*mm)) # size is (width, height)
for cond in randls[['iter' + str(x) for x in range(0,20)]]:
    sns.lineplot(ax=ax, data=randls, x='aa_pos', y=cond, zorder=2, alpha=1, legend=None)
ax.set(xlabel=None, xlim=(1,1616), xticks=[], ylabel=None, ylim=(-12,12), yticks=[])
ax.grid(b=None)
color_domains(ax=ax, domains=doms_dnmt, palette=pal_dnmt, spines=False)
fig.set_size_inches(32*mm, 24*mm) # sometimes need to resize due to monitor issues
# plt.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Fig2_schematic_1d_clus_randomized_loess_v1_211217.pdf', **plot_kws_savepdf)

### Fig 2a schematic -- observed LOESS lineplot w/ 95% CI of the null LOESS distribution as shaded region
fig, ax = plt.subplots(figsize=(32*mm, 24*mm)) # size is (width, height)
sns.lineplot(ax=ax, data=df_1d, x='aa_pos', y='func_56d_ls', zorder=2, **plot_kws_line)
ax.fill_between(x=df_1d['aa_pos'], y1=df_1d['95ci_min'], y2=df_1d['95ci_max'], alpha=0.4, color='#525252', ec='none', zorder=1.5)
ax.set(xlabel=None, xlim=(1,1616), xticks=[], ylabel=None, ylim=(-12,12), yticks=[])
ax.grid(b=None)
color_domains(ax=ax, domains=doms_dnmt, palette=pal_dnmt, spines=False)
ax.axvspan(xmin=544.5, xmax=545.5, color='#ffde5e', alpha=1)
fig.set_size_inches(32*mm, 24*mm) # sometimes need to resize due to monitor issues
# plt.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Fig2_schematic_1d_clus_obs_vs_null_loess_v1_211217.pdf', **plot_kws_savepdf)

### Fig 2a schematic -- histogram to compare obs. LOESS to null distribution of LOESS for single AA (aa 545)
fig, ax = plt.subplots(figsize=(32*mm, 24*mm)) # size is (width, height)
sns.kdeplot(ax=ax, data=randls[randls['aa_pos'] == 545][randls.columns[2:]].transpose(), palette=['#525252'], fill=True, legend=None,)# cumulative=True)
ax.set(xlabel='LOESS score', xticks=[], ylabel='Density', yticks=[])
ax.grid(b=None)
ax.axvline(x=randls[randls['aa_pos'] == 545]['func_56d_ls'].values[0], color='#e41a1c', alpha=1)
fig.set_size_inches(32*mm, 24*mm) # sometimes need to resize due to monitor issues
# plt.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Fig2_schematic_1d_clus_aa_loess_dist_histogram_v1_211219.pdf', **plot_kws_savepdf)


### Fig 2b -- line plot of -log10(corrected p-value) vs. CDS position
fig, ax = plt.subplots(figsize=(40*mm, 30*mm)) # size is (width, height)
sns.lineplot(ax=ax, data=df_1d, x='aa_pos', y='log10', zorder=1.5, color='#525252', alpha=0.8, legend=None)
# sns.scatterplot(ax=ax, data=df_1d[df_1d['sig']], x='aa_pos', y='log10', color='#e41a1c', alpha=1, s=9, ec='none', clip_on=False, legend=None, zorder=1.6)
ax.set(xlabel=None, xlim=(1,1616), xticks=[1,1616], ylabel='-log10(p-value)', ylim=(0,2), yticks=[0,1,2])
ax.grid(b=None)
ax.axhline(y=(-1 * np.log10(0.05)), c='k', alpha=0.7, ls='--')
color_domains(ax=ax, domains=doms_dnmt, palette=pal_dnmt, spines=False)
fig.set_size_inches(60*mm, 45*mm) # sometimes need to resize due to monitor issues
# plt.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Fig2_1d_clus_pval_cds_lineplot_v1_211220.pdf', **plot_kws_savepdf)


### Fig 2c -- clustered heatmap of gRNAs in PWES analysis
pal_clus = {k:'#CCCCCC' for k in range(1,20)}
pal_clus[1], pal_clus[2] = '#e41a1c', '#377eb8'
plot_clus_clrs = df_clus['cl_new_rank'].apply(lambda x: pal_clus[x]).copy()
plot_clus_clrs.index = df_pws1.index.copy()
plot_clus_clrs._name = 'Cluster'

g = sns.clustermap(data=df_pws1, row_linkage=link, col_linkage=link, figsize=(80*mm, 80*mm), vmin=-1, vmax=1, cmap='RdBu_r', row_colors=plot_clus_clrs, colors_ratio=0.02,
                   xticklabels=False, yticklabels=False, cbar_pos=(0.05,0.87,0.15,0.025), cbar_kws=dict(ticks=[-1,0,1], label='PWES', orientation='horizontal'))
g.ax_col_dendrogram.set_visible(False)
g.ax_heatmap.set(xlabel=None, ylabel=None, rasterized=True)
for spine in g.ax_heatmap.spines.items():
    spine[1].set_visible(True)
    spine[1].set_color('k')
for spine in g.ax_cbar.spines:
    g.ax_cbar.spines[spine].set_color('k')
    g.ax_cbar.spines[spine].set_linewidth(0.5)
    g.ax_cbar.spines[spine].set_visible(True)
# mark the cluster boundaries
plot_cum = 0
for i in range(1,20):
    plot_cum += df_clus['cl_new'].value_counts()[i]
    if i != 19: g.ax_heatmap.axhline(y=plot_cum - 1, c='k', lw=0.5)
        # g.ax_row_colors.axhline(y=plot_cum-1, c='k', lw=0.5)
g._figure.set_size_inches(80*mm,80*mm)
# g.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Fig2_CHUNKS_3d_cluster_heatmap_v2_211220.pdf', **plot_kws_savepdf, bbox_inches=None, dpi=1200)
# g.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Fig2_CHUNKS_3d_cluster_heatmap_v1_211220.png', **plot_kws_savepng, bbox_inches=None)

### Fig 2d,e -- scatter plot of gRNA CDS cut sites and box plot of gRNA func scores by cluster
plot_kws_clscatter = dict(s=6.25, ec='none', clip_on=False, legend=None)
plot_kws_box = dict(boxprops=dict(facecolor='none', lw=0.5, alpha=1), showfliers=False, orient='h', order=np.arange(1,20))
plot_kws_strip = dict(dodge=True, ec='none', jitter=0.15, size=2, alpha=0.2, orient='h', order=np.arange(1,20))
pal_clus = {k:v for k,v in zip(range(1,20), ['#e41a1c','#377eb8'] + ['#525252' for i in range(17)])}

fig, (ax1,ax2) = plt.subplots(ncols=2, gridspec_kw={'width_ratios': [3,1]}, figsize=(100*mm, 50*mm))
sns.scatterplot(data=df_clus[df_clus['cl_new_rank'].isin([1,2])], x='aa_pos', y='cl_new_rank', ax=ax1, hue='cl_new_rank', palette=pal_clus, alpha=0.8, **plot_kws_clscatter)
sns.scatterplot(data=df_clus[~df_clus['cl_new_rank'].isin([1,2])], x='aa_pos', y='cl_new_rank', ax=ax1, color='#525252', alpha=0.6, **plot_kws_clscatter)
sns.boxplot(data=df_clus, x='func_56d', y='cl_new_rank', ax=ax2, palette=pal_clus, width=0.6, zorder=2, **plot_kws_box)
sns.stripplot(data=df_clus, x='func_56d', y='cl_new_rank', ax=ax2, palette=pal_clus, zorder=1, **plot_kws_strip)
for (i, plot_artist), plot_clr in zip(enumerate(ax2.artists), pal_clus.values()):
    plot_artist.set_edgecolor(plot_clr)
    for j in range(i*5, i*5 + 5):
        ax2.lines[j].set_linewidth(0.5)
        ax2.lines[j].set_color(plot_clr)
        ax2.lines[j].set_alpha(1)
ax1.grid(b=None), ax2.grid(b=None)
ax1.set(xlabel=None, ylabel='Cluster', xlim=(351,1600), xticks=[351,1600], ylim=(0.5,19.5), yticks=[x for x in range (1,20)], yticklabels=[x for x in range (1,20)])
ax2.set(xlabel='Resistance score', xlim=(-12,12), xticks=[-12,-6,0,6,12], ylim=(0.5 - 1, 19.5 - 1), yticks=[x - 1 for x in range (1,20)], yticklabels=[], ylabel=None)
for ax in [ax1,ax2]:
    ax.spines['top'].set_visible(True), ax.spines['right'].set_visible(True)
color_domains(ax1, doms_dnmt, pal_dnmt, dom_kws=dict(alpha=0.2), nodom_kws=dict(fc='#CCCCCC', alpha=0.2))
ax1.axhline(y=2.5, ls='--', lw=0.5, c='k', alpha=0.8)
ax2.axhline(y=2.5 - 1, ls='--', lw=0.5, c='k', alpha=0.8)
ax1.add_patch(mpl.patches.Rectangle(xy=(351,-1.5), width=1600-351, height=0.6, fc='#CCCCCC', ec='k', clip_on=False, zorder=1, alpha=0.8, lw=0.5))
ax1.add_patch(mpl.patches.Rectangle(xy=(518,-1.5), width=54, height=0.6, fc='#ffde5e', ec='k', clip_on=False, zorder=1.2, lw=0.5))
ax1.add_patch(mpl.patches.Rectangle(xy=(652,-1.5), width=50, height=0.6, fc='#ffde5e', ec='k', clip_on=False, zorder=1.2, lw=0.5))
fig.set_size_inches(100*mm, 50*mm)
# plt.savefig('./KN10164_clustering_analysis_v1_211201/plots/' + 'Fig2_3d_clusters_CDS_scatter_boxplot_v2_211221.pdf', **plot_kws_savepdf)

#%% 12DEC21: perform CHUNKS on gRNAs/AAs w/ revised conds (don't run again; import results)

### for gRNA clustering, find all pairwise gRNA sums and z-score (incl. non-resolved gRNAs)
### then tanh-scale and scale w/ gaussian (std=16), call clusters w/ cophen_dist t=13.8
### 19 clusters output -- renumber by mean func score rank

### calculate PWES score matrix for gRNAs
# calculate pw summed LFCs for all gRNAs targeting resolved residues in 4wxx (use standard rounding)
df_scores = df_dnmt.copy()
df_scores['aa_pos'] = df_scores['cut_site_AA'].round() # round grna cut sites
list_aas = df_scores[df_scores['aa_pos'].isin(df_gauss.index)]['aa_pos'] # get all resolved gRNAs (n=646)
# get pairwise sums of gRNA scores for resolved gRNAs
df_pws = df_scores[df_scores['aa_pos'].isin(list_aas)].copy() # get resolved gRNAs (n=646)
df_pws = pd.DataFrame(index=df_pws['sgRNA_ID'], columns=df_pws['sgRNA_ID'], data=df_pws['func_56d'].values[:, None] + df_pws['func_56d'].values[None, :]) # pairwise sums
# get pairwise sums of gRNA scores for ALL gRNAs (n=830) to find mean and stdev for z-scoring and tanh-scaling
temp = pd.DataFrame(index=df_scores['sgRNA_ID'], columns=df_scores['sgRNA_ID'], data=df_scores['func_56d'].values[:, None] + df_scores['func_56d'].values[None, :])
temp = pd.DataFrame(index=temp.index, columns=temp.columns, data=np.where(np.triu(np.ones(temp.shape), k=1).astype(bool), temp, np.nan))
temp2 = pd.Series([y for x in temp.columns for y in temp[x]], name='sum_lfc').dropna() # pw summed gRNA func scores
# z-score pairwise sums w/ pairwise sum mean and std of all gRNAs, then tanh-scale
df_pws1 = np.tanh((df_pws - temp2.mean()) / temp2.std()) # all gRNAs (n=830) pw_sum mean == -0.845; std == 4.406
df_pws1.index, df_pws1.columns = list_aas, list_aas # replace index/columns
df_pws1 = df_pws1 * df_gauss.loc[list_aas, list_aas].copy()

### perform hierarchical clustering for gRNAs PWES scores
# import old clustering data for comparison
df_clus = df_scores.loc[df_scores['aa_pos'].isin(list_aas)][list_refcols + ['func_56d','enrich_56d','aa_pos']].copy()
df_clus['bl106'] = np.where(df_clus['sgRNA_ID'].isin(list_bl106), True, False)
df_clus = df_clus.merge(pd.read_csv('./KN10164_clustering_analysis_v1_211201/KN10152_DNMT1_clusters_v1_200722.csv', index_col=0)[['sgRNA_ID','Cluster']],
                        how='left', on='sgRNA_ID').rename(columns={'Cluster':'cl_old'})

# cluster using cophen_dist(t)=13.8 to cluster (19 clusters output)
link = sp_cl.hierarchy.linkage(df_pws1, method='ward', metric='euclidean', optimal_ordering=True)
df_clus['cl_new'] = sp_cl.hierarchy.fcluster(link, t=13.8, criterion='distance')
# re-organize clusters by func_56d mean of their component gRNAs (rank 1 == old cluster 1, rank 2 == old cluster 10 (partial))
list_cl_ranks = df_clus.groupby('cl_new')['func_56d'].mean().rank(ascending=False)
df_clus['cl_new_rank'] = df_clus['cl_new'].apply(lambda x: list_cl_ranks[x])
# get cluster stats (func_56d mean and # of gRNAs per cluster)
df_clstats = df_clus.groupby(by='cl_old', as_index=False).agg({'func_56d':'mean', 'cl_old':'count'}).rename(columns={'func_56d':'cl_old_mean'})
df_clstats = df_clstats.merge(df_clus.groupby(by='cl_new', as_index=False).agg({'func_56d':'mean', 'cl_new':'count'}), how='outer', left_index=True, right_index=True).rename(columns={'func_56d':'cl_new_mean'})
df_clstats = df_clstats.merge(df_clus.groupby(by='cl_new_rank', as_index=False).agg({'func_56d':'mean', 'cl_new_rank':'count'}), how='outer', left_index=True, right_index=True).rename(columns={'func_56d':'cl_rank_mean'})
df_clstats['cluster'] = range(1, len(df_clstats.index) + 1)
df_clstats = df_clstats.reindex(columns=['cluster','cl_old','cl_old_mean','cl_new','cl_new_mean','cl_new_rank','cl_rank_mean'])

### calculate empirical p-vals for the mean gRNA scores per cluster
temp_pv = {}
for cl, cl_idxs in df_clus.groupby('cl_new_rank').groups.items():
    cl_grnas = df_clus.loc[cl_idxs, 'sgRNA_ID'].copy()
    temp_dist = df_rand.loc[df_rand['sgRNA_ID'].isin(cl_grnas), df_rand.columns[3:]].mean()
    temp_pv[cl] = (10000 - temp_dist.lt(df_clstats[df_clstats['cluster'] == cl]['cl_rank_mean'].values[0]).sum()) / 10000
df_clstats = df_clstats.merge(pd.DataFrame.from_dict(orient='index', columns=['cl_rank_pval'], data=temp_pv), how='outer', left_on='cluster', right_index=True)
# BH FDR correction of p-values
temp_corr = smm.multipletests(df_clstats['cl_rank_pval'], alpha=0.05, method='fdr_bh') # apply benjamini-hochberg FDR correction
df_clstats['sig'] = temp_corr[0]
df_clstats['corr_pval'] = temp_corr[1]

### export the gRNAs 3D clustering results and stats to csv
### df_clus.to_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_grna_clusters_v1_211213.csv', index=False)
### df_clstats.to_csv('./KN10164_clustering_analysis_v1_211201/' + 'dnmt1_chunks_grna_clusters_stats_v1_211213.csv', index=False)

