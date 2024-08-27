

# calculate distances

# calculate pairwise scores

# combine and hierarchical clustering


from pathlib import Path

import numpy as np
import pandas as pd

import scipy.cluster as sp_cl
import statsmodels.stats.multitest as smm

from be_scan.plot.clustering_helper import *


# defining lists of common column headers for downstream ease
list_refcols = ['sgRNA_ID', 'sgRNA_seq', 'Gene', 'cut_site_NT', 'cut_site_AA', 'Domain']
list_condcols = ['DAC_res', 'GSKi_res']
list_enrichcols = ['enrich_DAC', 'enrich_GSKi']

# subdirectory to hold analysis/output files (05OCT21)
outpath = Path.cwd() /'KN10164_analysis_v1_211004'

### importing processed CRISPR screen data for DNMT1/UHRF1 ###
in_lfc = pd.read_csv('./DNMT1_screen_processing/DNMT1_screen_master_v1.csv')
df_lfc = pd.read_csv('./DNMT1_screen_processing/DNMT1_screen_reduced_v1.csv')
df_dnmt = pd.read_csv('./DNMT1_screen_processing/DNMT1_screen_reduced_dnmt1_v1.csv',)
df_uhrf = pd.read_csv('./DNMT1_screen_processing/DNMT1_screen_reduced_uhrf1_v1.csv')

# identify enriched sgRNAs across both screens
df_enrich = in_lfc[(in_lfc['enrich_DAC'] == 'enrich') | (in_lfc['enrich_GSKi'] == 'enrich')].copy()

#%% 1D clustering analysis
### perform LOESS on DNMT1 CRISPR screen data w/ params: span=100 AA/L, interp=quadratic

df_dnmtls = pd.DataFrame(columns=['aa_pos'])
df_dnmtls['aa_pos'] = np.arange(1,1617)
df_dnmtls['DAC_res_ls'] = loess_v3(x_obs=df_dnmt['cut_site_AA'], y_obs=df_dnmt['DAC_res'],
                                   x_out=df_dnmtls['aa_pos'], span=100/1616, interp_how='quadratic')['y_loess']

# calculate null distribution for each AA by shuffling gRNA function scores and running LOESS
# only run this section once (takes forever); afterwards, just import exported files

df_rand = df_dnmt[['sgRNA_ID','cut_site_AA','DAC_res']].copy()
df_randls = df_dnmtls[['aa_pos','DAC_res_ls']].copy()
for i in range(10000):
    df_rand['iter' + str(i)] = df_rand['DAC_res'].sample(frac=1, ignore_index=True)
    df_randls['iter' + str(i)] = loess_v3(x_obs=df_rand['cut_site_AA'], y_obs=df_rand['iter' + str(i)],
                                          x_out=df_randls['aa_pos'], span=100/1616, interp_how='quadratic')['y_loess']

### export null distributions for shuffled gRNA scores and AA loess scores for future import bc it takes forever
# df_rand.to_csv('./DNMT1_clustering_analysis/dnmt1_grna_DAC_res_null_dist.csv', index=False)
# df_randls.to_csv('./DNMT1_clustering_analysis/dnmt1_aa_loess_DAC_res_null_dist.csv', index=False)

# calculate mean, std, 95% CI, etc. for the AA LOESS values
df_1d = df_randls[['aa_pos','DAC_res_ls']].copy()
df_1d['mean'] = df_randls.iloc[:,2:].mean(axis=1)
df_1d['std'] = df_randls.iloc[:,2:].std(axis=1)
df_1d['95ci_min'] = df_randls.iloc[:,2:].quantile(q=0.025, axis=1)
df_1d['95ci_max'] = df_randls.iloc[:,2:].quantile(q=0.975, axis=1)
df_1d['obs_gt'] = df_randls.iloc[:,2:].gt(df_randls['DAC_res_ls'], axis=0).sum(axis=1) # get # of values greater than obs_val
df_1d['1t_pval'] = df_1d['obs_gt'] / 10000 # divide "rank" of obs val by N to get empirical p-val

temp = smm.multipletests(df_1d['1t_pval'], alpha=0.05, method='fdr_bh') # apply benjamini-hochberg FDR correction
df_1d['sig'] = temp[0]
df_1d['corr_pval'] = temp[1]
df_1d['log10'] = np.log10(df_1d['corr_pval']) * -1 # get -log10(p-values) for plotting

### export the obs vs. null statistics and pvals
# df_1d.to_csv('./DNMT1_clustering_analysis/dnmt1_1d_clustering_stats.csv', index=False)

#%% 3D clustering analysis

### get centroids and calculate pairwise distances (only run once, then import)
df_cent = get_centroids(pdb_id='4wxx', aa_sel='4wxx and chain B and polymer.protein',
                        save_centroids=True, out_csv='dnmt1_4wxx_centroids.csv')
df_pwdist = get_pairwise_dist(df_cent)
df_pwdist.to_csv('./DNMT1_clustering_analysis/4wxx_pairwise_dist.csv')

### importing centroids/pairwise distances
df_cent = pd.read_csv('./DNMT1_clustering_analysis/dnmt1_4wxx_centroids.csv')
df_pwdist = pd.read_csv('./DNMT1_clustering_analysis/4wxx_pairwise_dist.csv', index_col=0)
df_pwdist.columns = df_pwdist.columns.astype(int) # make columns int64 dtype

### scale pairwise distances with gaussian fx (std=16 angstroms)
df_gauss = df_pwdist.apply(lambda x: gauss(x, std=16)).copy()

### calculating pairwise summed sgRNA resistance scores
# make pairwise sums matrix for all DNMT1 sgRNAs and resolved DNMT1 sgRNAs
df_scores = df_dnmt.copy()
df_scores['aa_pos'] = df_scores['cut_site_AA'].round() # round grna cut sites
list_aas = df_scores[df_scores['aa_pos'].isin(df_gauss.index)]['aa_pos'] # get all resolved gRNAs (n=646)
df_pws = df_scores[df_scores['aa_pos'].isin(list_aas)].copy() # pairwise matrix for resolved sgRNAs
df_pws = pd.DataFrame(index=df_pws['sgRNA_ID'], columns=df_pws['sgRNA_ID'],
                      data=df_pws['DAC_res'].values[:, None] + df_pws['DAC_res'].values[None, :]) # pairwise sums matrix
# calculate pairwise sums for ALL DNMT1 sgRNAs (n=830) to find mean and stdev for z-scoring and tanh-scaling
df_pws_all = pd.DataFrame(index=df_scores['sgRNA_ID'], columns=df_scores['sgRNA_ID'],
                           data=df_scores['DAC_res'].values[:, None] + df_scores['DAC_res'].values[None, :])
df_pws_all = pd.DataFrame(index=df_pws_all.index, columns=df_pws_all.columns,
                          data=np.where(np.triu(np.ones(df_pws_all.shape), k=1).astype(bool), df_pws_all, np.nan))
# collapse pw summed scores into series to calculate mean and std for z-scoring/tanh-scaling
list_pws_all = pd.Series([y for x in df_pws_all.columns for y in df_pws_all[x]], name='sum_lfc').dropna()
# z-score resolved sgRNA scores (df_pws) w/ pairwise sum mean and std of all sgRNAs, then tanh-scale the z-scores
df_pws = np.tanh((df_pws - list_pws_all.mean()) / list_pws_all.std()) # all gRNAs (n=830) pw_sum mean == -0.845; std == 4.406
df_pws.index, df_pws.columns = list_aas, list_aas # replace index/columns with amino acid positions
# scale the tanh-normalized pairwise scores by the distance component (gaussian)
df_pws = df_pws * df_gauss.loc[list_aas, list_aas].copy()

# perform hierarchical clustering using cophen_dist(t)=13.8 to cluster (19 clusters output)
link = sp_cl.hierarchy.linkage(df_pws, method='ward', metric='euclidean', optimal_ordering=True)
df_clusters = df_scores.loc[df_scores['aa_pos'].isin(list_aas)][list_refcols + ['aa_pos','DAC_res','enrich_DAC']].copy()
df_clusters['clusters_out'] = sp_cl.hierarchy.fcluster(link, t=13.8, criterion='distance')
# re-number clusters by mean DAC_res score of their component gRNAs (high to low)
list_cl_ranks = df_clusters.groupby('clusters_out')['DAC_res'].mean().rank(ascending=False)
df_clusters['clusters_rank'] = df_clusters['clusters_out'].apply(lambda x: list_cl_ranks[x])

# get cluster stats (DAC_res mean and # of gRNAs per cluster)
df_cl_stats = df_clusters.groupby(by='clusters_rank', as_index=False).agg({'DAC_res':'mean', 'clusters_out':'count'}).rename(columns={'DAC_res':'cluster_mean'})
df_cl_stats['cluster_id'] = range(1, len(df_cl_stats.index) + 1)
df_cl_stats = df_cl_stats.reindex(columns=['cluster_id','cluster_mean','clusters_rank','clusters_out'])

### calculate empirical p-vals for the mean gRNA scores per cluster
# reusing the shuffled sgRNA resistance scores null dist from earlier
temp_pvals = {}
for cluster, cluster_idxs in df_clusters.groupby('clusters_rank').groups.items():
    cluster_grnas = df_clusters.loc[cluster_idxs, 'sgRNA_ID'].copy()
    temp_dist = df_rand.loc[df_rand['sgRNA_ID'].isin(cluster_grnas), df_rand.columns[3:]].mean()
    temp_pvals[cluster] = (10000 - temp_dist.lt(df_cl_stats[df_cl_stats['cluster_id'] == cluster]['cluster_mean'].values[0]).sum()) / 10000
df_cl_stats = df_cl_stats.merge(pd.DataFrame.from_dict(orient='index', columns=['cluster_pval'], data=temp_pvals), how='outer', left_on='cluster_id', right_index=True)
# BH FDR correction of p-values
temp_corr = smm.multipletests(df_cl_stats['cluster_pval'], alpha=0.05, method='fdr_bh')
df_cl_stats['sig'] = temp_corr[0]
df_cl_stats['corr_pval'] = temp_corr[1]

#%% Supp. Fig. 2 -- weighted average proximity (WAP) scores and cluster significance
### WAP = sum of all pairwise enrichment scores (absolute value) in a given cluster
### simulate null dist of WAPs by shuffling AA pos targeted by resolved sgRNAs
### and recalculating the PWES and WAP scores

# generate shuffled aa_positions for resolved sgRNAs (n=10000)
# only run once and import afterwards
df_sim = df_scores[df_scores['aa_pos'].isin(df_gauss.index)][['sgRNA_ID','DAC_res','aa_pos']].copy() # get all resolved sgRNAs (n=646)
df_sim = df_sim.rename(columns={'aa_pos':'aa_obs'}).reset_index(drop=True) 
temp_dict = {}
for i in np.arange(0,10000):
    temp_dict['iter' + str(i)] = df_sim['aa_obs'].sample(frac=1, ignore_index=True).rename('iter' + str(i))
df_sim = pd.concat([df_sim] + [x for x in temp_dict.values()], axis=1)
# export to csv for future import
# df_sim.to_csv('./DNMT1_clustering_analysis/dnmt1_3d_clustering_shuffled_aa_pos.csv', index=False)

### calculate the z-scored/tanh-transformed pairwise sums matrix
# make pairwise sums matrix for resolved sgRNAs
temp_pws = pd.DataFrame(index=df_sim['sgRNA_ID'], columns=df_sim['sgRNA_ID'],
                        data=df_sim['DAC_res'].values[:, None] + df_sim['DAC_res'].values[None, :])
# calculate pairwise sums for ALL DNMT1 sgRNAs (n=830) to find mean and stdev for z-scoring and tanh-scaling
temp = pd.DataFrame(index=df_scores['sgRNA_ID'], columns=df_scores['sgRNA_ID'], data=df_scores['DAC_res'].values[:, None] + df_scores['DAC_res'].values[None, :])
temp = pd.DataFrame(index=temp.index, columns=temp.columns, data=np.where(np.triu(np.ones(temp.shape), k=1).astype(bool), temp, np.nan))
temp2 = pd.Series([y for x in temp.columns for y in temp[x]], name='sum_lfc').dropna() # pw summed gRNA func scores
# z-score resolved sgRNA scores (df_pws) w/ pairwise sum mean and std of all sgRNAs, then tanh-scale the z-scores
temp_pws = np.tanh((temp_pws - temp2.mean()) / temp2.std()) # all gRNAs (n=830) pw_sum mean == -0.845; std == 4.406
### create a 646x646 mask to set self vs. self diagonal to 0
temp_mask = np.tril(np.ones(temp_pws.shape)).astype(bool)
### make dict of (cluster, list of sgRNAs)
dict_clus = {k: df_clusters[df_clusters['cluster_id'] == k]['sgRNA_ID'].tolist() for k in np.arange(1,20)}


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
# export simulated WAP null dist for future import
# df_wap.to_csv('./DNMT1_clustering_analysis/dnmt1_3d_clustering_shuffled_aa_WAP_null_dist.csv', index=False)

### calculate stats/p-vals for observed WAP vs. null distribution
df_wap_stats = df_wap[['cluster','obs_wap']].copy()
df_wap_stats['mean'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].mean(axis=1)
df_wap_stats['std'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].std(axis=1)
df_wap_stats['95ci_min'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].quantile(q=0.025, axis=1)
df_wap_stats['95ci_max'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].quantile(q=0.975, axis=1)
df_wap_stats['obs_gt'] = df_wap[['iter' + str(i) for i in np.arange(0,10000)]].gt(df_wap['obs_wap'], axis=0).sum(axis=1)
df_wap_stats['1t_pval'] = df_wap_stats['obs_gt'] / 10000
# df_wap_stats.to_csv('./DNMT1_clustering_analysis/dnmt1_3d_clustering_WAP_stats.csv', index=False)
