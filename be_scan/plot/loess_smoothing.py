"""
Author: Calvin XiaoYang Hu
Adapted from: Ceejay Lee - 1D_LOESS_clustering.ipynb
Date: 240313

{Description: This set of functions calculates the loess smoothed values across thelength of the gene, 
              calculates a randomized rearrangement as a baseline, 
              and then determines the significance of smoothed hits}
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import random

# specific to the loess function
import scipy.interpolate as interp
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.stats.multitest as smm
from matplotlib.collections import LineCollection

interp_methods = ['linear', 'nearest', 'nearest-up', 'zero', 
                  'slinear', 'quadratic', 'cubic', 'previous', 'next']

def loess_smoothing(
    df, x_column, comparisons, span, 
    interp_method='quadratic', n_repeats=1000,
    loess_kws=None, interp_kws=None, smm_multipletests_kws=None,
):
    assert interp_method in interp_methods

    # OR OPERATOR ASSIGNS VALUE IF NOT NONE, ELSE ASSIGNS A DEFAULT #
    loess_kws = loess_kws or {'missing':'raise', 'return_sorted':False, 'it':0}
    interp_kws = interp_kws or {'fill_value':'extrapolate'}
    smm_multipletests_kws = smm_multipletests_kws or {
        'alpha': 0.05, 'method': 'fdr_bh', 'is_sorted': False, 'returnsorted': False
    }

    # FILTER AND FLOOR DATA #
    df_filtered = df[[x_column] + comparisons].copy()
    df_filtered = df_filtered[df_filtered[x_column] >= 0].sort_values(by=x_column)
    df_filtered[comparisons] = df_filtered[comparisons].apply(np.floor)

    results = {}
    for comp in comparisons:
        x_obs, y_obs = df_filtered[x_column], df_filtered[comp]
        x_min, x_max = x_obs.min(), x_obs.max()
        span_new = span / (x_max - x_min) if span > 1 else span

        df_loess = loess_v3(x_obs, y_obs, span=span_new,
                            interp_method=interp_method,
                            loess_kws=loess_kws, interp_kws=interp_kws)
        df_rand = randomize(x_obs, y_obs, span=span_new,
                            n_repeats=n_repeats,
                            interp_method=interp_method, random_seed=0)
        df_pvals = calculate_sig(df_loess, df_rand,
                                 n_repeats=n_repeats,
                                 smm_multipletests_kws=smm_multipletests_kws)

        results[comp] = {
            'loess': df_loess, 'rand': df_rand, 'pvals': df_pvals
        }
    return results

def plot_loess(
    results, gene, 
    domains=[], domains_alpha=0.25, 
    show=True, savefig=False, out_name='loesssmoothing', out_type='png', out_dir='',
):
    for comp, res in results.items():
        df_pvals = res['pvals']
        df_loess = res['loess']

        # --- P-VALUE PLOT --- #
        fig1, ax1 = plt.subplots(figsize=(6, 4))
        df_plotting = df_pvals[['corr_pval']].copy()
        df_plotting['-log10'] = np.log10(df_plotting['corr_pval'] + 10**-4) * -1
        # X Y VALUES #
        x, y = df_plotting.index.to_numpy(), df_plotting['-log10'].to_numpy()
        x, y = increase_resolution(x), increase_resolution(y)
        threshold = np.log10(0.05 + 10**-4) * -1
        # PLOT ABOVE AND BELOW DIFFERENT COLORS #
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        colors = ['red' if (yi > threshold and yj > threshold) else 'lightgray' for yi, yj in zip(y[:-1], y[1:])]
        line = LineCollection(segments, colors=colors, linewidths=1)
        ax1.add_collection(line)
        # PLOT SETTINGS #
        ax1.axhline(y=threshold, ls='--', c='k', linewidth=1)
        ax1.set_title(f'{gene} Linear Clustering')
        ax1.set_ylabel('-log 10 (adjusted P value)')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['bottom'].set_visible(False)
        plt.xticks([])
        # DOMAINS #
        for d in domains:
            ax1.axvspan(d['start'], d['end'], alpha=domains_alpha, facecolor=d['color'])
        # SAVE AND CLOSE #
        if savefig:
            plt.savefig(Path(out_dir) / f'{out_name}_{comp}_pvals.{out_type}', dpi=300)
        if show: plt.show()
        plt.close()

        # --- LOESS SIGNAL PLOT --- #
        fig2, ax2 = plt.subplots(figsize=(6, 4))
        df_plotting = df_loess[['y_loess']].copy()
        # PLOT #
        ax2.plot(df_plotting.index, df_plotting['y_loess'], color='gray', markersize=2)
        # PLOT SETTINGS #
        ax2.set_title(f'{gene} Loess Smoothed')
        ax2.set_ylabel('sgRNA Score')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.axhline(y=0, ls='--', c='k', linewidth=1)
        plt.xticks([])
        # DOMAINS #
        for d in domains:
            ax2.axvspan(d['start'], d['end'], alpha=domains_alpha, facecolor=d['color'])
        # SAVE AND CLOSE #
        if savefig:
            plt.savefig(Path(out_dir) / f'{out_name}_{comp}_loess.{out_type}', dpi=300)
        if show: plt.show()
        plt.close()

# HELPER FUNCTIONS #

def loess_v3(x_obs, y_obs, span, interp_method, 
    loess_kws={}, interp_kws={}, 
    ):

    """[Summary]
    this function calculates the smoothed arbitrary values for each value 
    in your x_out for the enrichment profile of your POI

    Returns
    ----------
    df_loess: pandas dataframe
        dataframe of x_vals and y_obs from input, 
        and y_loess and type from output corresponding to 
        the loess calculated or interpolated values of the loess smoothing and the type
    """
    x_obs = x_obs.astype(float)
    # make the resolution the same as x values
    df_loess = pd.DataFrame()
    df_loess['x_vals'] = np.arange(x_obs.min(), x_obs.max()+1, 1)

    loess_out = lowess(endog=y_obs, exog=x_obs, frac=span, **loess_kws)
    if interp_method == 'statsmodels': 
        df_loess['y_loess'] = loess_out
    else: 
        # for some reason sometimes loess_out is 2D and sometimes its 1D
        loess_out = np.array(loess_out)
        df_interp = pd.DataFrame(loess_out)
        if df_interp.shape[1] == 2: 
            df_interp = df_interp.drop(columns=[0])
            df_interp = df_interp.rename({1:0}, axis='columns')
        df_interp['x_obs'] = x_obs.values
        df_interp = df_interp.rename({0: 'y_obs_loess'}, axis='columns')
        # in case theres duplicate values in x_obs
        df_interp = df_interp.groupby('x_obs', as_index=False).agg({'y_obs_loess':'mean'}) 
        # interpolate missing variables
        fx_interp = interp.interp1d(x=df_interp['x_obs'], y=df_interp['y_obs_loess'], 
                                    kind=interp_method, **interp_kws)
        # apply loess interpolated function to a range of x values
        df_loess['y_loess'] = fx_interp(df_loess['x_vals'])

    df_loess['type'] = np.where(df_loess['x_vals'].isin(x_obs), 'loess', 'interp')
    return df_loess

def randomize(x_obs, y_obs, span, 
    n_repeats=1000, interp_method='quadratic', 
    ): 

    """[Summary]
    this function provides the null distribution to compare your actual loess again; 
    it calculates the smoothed arbitrary values for each value 
    in your x_out (see below) for the shuffled enrichment profile of your POI

    Returns
    ----------
    df_rand: pandas dataframe
        has n_repeats number of columns each of calculated loess scores from randomized y_obs
    """
    loess_args = {'x_obs':x_obs, 'y_obs':random.sample(list(y_obs), len(y_obs)), 
                  'span':span, 'interp_method':interp_method}
    loess_out = loess_v3(**loess_args)

    loess_rand = []
    for _ in range(0, n_repeats): 
        loess_args['y_obs'] = random.sample(list(y_obs), len(y_obs))
        loess_out = loess_v3(**loess_args)
        loess_rand.append(loess_out.y_loess)

    df_rand = pd.concat([loess_out.x_vals] + loess_rand + [loess_out.type], 
                        ignore_index=True, axis=1 )
    return df_rand

def calculate_sig(df_loess, df_rand, 
    n_repeats=1000, smm_multipletests_kws={}
    ):

    """[Summary]
    this function calculates the significance of a "spike" in your loess-ed enrichment profile
    by calculating if it's > 95% of null distribution

    Returns
    ----------
    df_pvals: pandas dataframe
        dataframe with your position values (x_out) and corresponding loess values (y_loess) from input, 
        and empirical p-val (sig) as well as the corrected p-values (corr_pval)
    """

    df_pvals = pd.DataFrame({'x_vals': list(df_loess.x_vals), 
                             'y_loess': list(df_loess.y_loess)}, 
                            columns=['x_vals', 'y_loess'])
    # get num of values greater than y_loess
    # df_pvals['obs_gt'] = df_rand[df_rand.columns[1:-1]].gt(df_pvals['y_loess'], axis=0).sum(axis=1)
    gt_counts = np.sum(df_rand[df_rand.columns[1:-1]] > df_loess['y_loess'].values[:, None], axis=1)
    # divide "rank" of obs val by N to get empirical p-val
    df_pvals['1t_pval'] = gt_counts / n_repeats 

    # apply benjamini-hochberg FDR correction
    temp = smm.multipletests(df_pvals['1t_pval'], **smm_multipletests_kws)
    df_pvals['sig'] = temp[0]
    df_pvals['corr_pval'] = temp[1]
    return df_pvals

def output_clusters(name, xvals, yvals): 
    clusters_x = []
    for x, y in zip(xvals, yvals): 
        if y > 1.30103: 
            clusters_x.append(x)
    if len(clusters_x) > 0: 
        diffs = np.diff(clusters_x)
        discontinuities = np.where(diffs > 1)[0]
        breaks = np.concatenate(([0], discontinuities + 1, [len(clusters_x)]))
        ranges = [f"{clusters_x[start]}-{clusters_x[end-1]}" if start != end-1 else f"{clusters_x[start]}" 
                for start, end in zip(breaks[:-1], breaks[1:])]
        print(name, ':', ranges)

def increase_resolution(x, times=5): 
    for _ in range(times): 
        new_x = []
        for i in range(len(x) - 1):
            new_x.append(x[i])
            avg = (x[i] + x[i+1]) / 2
            new_x.append(avg)
        new_x.append(x[-1])
        x = new_x
    return x

# loess_smoothing(
#     df_filepath='tests/test_data/plot/NZL10196_v9_comparisons.csv', 
#     x_column='Edit_site_3A1', 
#     comparisons=["d3-pos", "d3-neg", "d6-pos", ], 
#     span=20, n_repeats=3000, 
#     domains=[{'start':300, 'end':350}], 
# )

# loess_range(
#     df_filepath='tests/test_data/plot/NZL10196_v9_comparisons.csv', 
#     x_column='Edit_site_3A1', 
#     comparison="d3-pos", 
#     spans=[10, 20, 30, 40], n_repeats=3000, 
#     domains=[{'start':300, 'end':350}], 
# )
