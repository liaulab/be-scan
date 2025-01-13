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
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path
import random

# specific to the loess function
import scipy.interpolate as interp
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.stats.multitest as smm

def loess_smoothing(df_filepath, 
                    x_column, comparisons, span, 

    interp_method='quadratic', n_repeats=1000, 
    savefig=True, show=True, out_name='loesssmoothing', out_type='png', out_dir='', return_df=True, # output params
    domains=[], domains_alpha=0.25, domains_color='lightblue', # draw domains

    subplots_kws={}, 
    loess_kws={'missing':'raise', 'return_sorted':False, 'it':0}, 
    interp_kws={'fill_value':'extrapolate'}, 
    smm_multipletests_kws={'alpha':0.05, 'method':'fdr_bh', 'is_sorted':False, 'returnsorted':False}, 
    ): 

    """[Summary]
    This function calculates the smoothed arbitrary values for each value 
    in your x_out for the enrichment profile of your POI
    
    Parameters
    ------------
    df_filepath : str, required
        filepath to .csv data
    x_column : str, required
        column of .csv, typically amino acid position
    comparisons : list of str, required
        list of comparisons that correspond to columns of data
    span : float, required (recommended range to start is ~0.05-0.10)
        hyperparameter for the span that you want to smooth over, requires optimization

    interp_method: str, optional, defaults to 'quadratic'
        method should 'statsmodels' or on be one of 
        'linear', 'nearest', 'nearest-up', 'zero', 
        'slinear', 'quadratic', 'cubic', 'previous', or 'next'
        from https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
    n_repeats: int, optional, defaults to 10,000
        the number of times you want to shuffle
        1000x is recommended for final data

    savefig : bool, optional, defaults to True
        whether or not to save the figure
    show : bool, optional, defaults to True
        whether or not to show the plot
    out_name : str, optional, defaults to 'loess_smoothed'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
        path to output directory
    return_df: bool, optional, defaults to True
        whether or not to return the dataframes

    domains : dict, optional, defaults to {}
        a list of start dicts of start and end for each domain annotation ie [{'start': 1, 'end': 2}, ]
    domains_alpha : float, optional, defaults to 0.25
        The level of transparency for the domains
    domains_color : str, optional, defaults to 'lightblue'
        The color for the domains

    subplots_kws : dict, optional, defaults to {}
`       input params for plt.subplots
        https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
    loess_kws: dict, optional, defaults to
        {'missing':'raise', 'return_sorted':False, 'it':0}
        input params for statsmodels.nonparametric.smoothers_lowess.lowess()
        https://www.statsmodels.org/devel/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html
    interp_kws: dict, optional, defaults to
        {'fill_value':'extrapolate'}
        input params for scipy.interpolate.interp1d()
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
    smm_multipletests_kws: dict, optional, defaults to 
        {'alpha':0.05, 'method':'fdr_bh', 'is_sorted':False, 'returnsorted':False}
        input params for smm.multipletests()
        https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html

    Returns
    ------------
    """
    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)
    
    # check variable inputs
    assert interp_method in ['linear', 'nearest', 'nearest-up', 'zero', 
                             'slinear', 'quadratic', 'cubic', 'previous', 'next']

    # check all column names are valid, and filter down data
    assert x_column in df_data.columns, "check param: x_column"
    for comp in comparisons: 
        assert comp in df_data.columns, f"check param: comparisons {comp}"

    df_filtered = df_data[[x_column] + comparisons]
    df_filtered = df_filtered[df_filtered[x_column] >= 0].sort_values(by=x_column)
    result = {}

    mpl.rcParams.update({'font.size': 10})
    fig, axes = plt.subplots(nrows=len(comparisons), ncols=1, 
                             figsize=(10, 2.5*len(comparisons)), **subplots_kws)
    if len(comparisons) == 1: axes = [axes]
    # process columns from dataframe and input into loess_v3
    for ax, comp in zip(axes, comparisons): 
        x_obs, y_obs = df_filtered[x_column], df_filtered[comp]

        # ADJUSTING SPAN #
        x_min, x_max = df_filtered[x_column].min(), df_filtered[x_column].max()
        if span > 1: 
            assert span < (x_max-x_min), 'Please enter a span less than the length of your protein.'
            span = span / (x_max-x_min) ### assumming continuous scanning

        # LOESS SMOOTHES ACROSS LENGTH OF PROTEIN FOR SIGNAL #
        df_loess = loess_v3(x_obs=x_obs, y_obs=y_obs, span=span, 
                            interp_method=interp_method, loess_kws=loess_kws, interp_kws=interp_kws, )
        # RANDOMIZES TO OBTAIN A BACKGROUND SIGNAL #
        df_rand = randomize(x_obs=x_obs, y_obs=y_obs, span=span, 
                            n_repeats=n_repeats, interp_method=interp_method, )
        # COMPARES SIGNAL AGAINST BACKGROUND SIGNAL #
        df_pvals = calculate_sig(df_loess=df_loess, df_rand=df_rand, 
                                 n_repeats=n_repeats, smm_multipletests_kws=smm_multipletests_kws, )

        df_plotting = df_pvals[['corr_pval']].copy()
        df_plotting['-log10'] = np.log10(df_plotting['corr_pval'] + 10**-4) * -1
        ax.plot(df_plotting.index, df_plotting['-log10'], color='steelblue', markersize=3)
        ax.axhline(y=np.log10(0.05 + 10**-4)*-1,ls='--', c='k', linewidth=1)
        ax.set_title(comp)
        for d in domains:
            ax.axvspan(d['start'], d['end'], alpha=domains_alpha, facecolor=domains_color)
        output_clusters(comp, df_plotting.index, df_plotting['-log10'])

        if return_df: result[comp] = {'loess':df_loess, 'rand':df_rand, 'pvals':df_pvals}

    plt.tight_layout() # ADJUST DIMENSIONS #
    # Save file
    outpath = Path(out_dir)
    if savefig: 
        out = f'{out_name}.{out_type}'
        plt.savefig(outpath / out, format=out_type, dpi=300)
    if show: plt.show()
    plt.close()

    if return_df: return result
    else: return None


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

# loess_smoothing(
#     df_filepath='tests/test_data/plot/NZL10196_v9_comparisons.csv', 
#     x_column='Edit_site_3A1', 
#     comparisons=["d3-pos", "d3-neg", "d6-pos", ], 
#     span=0.05, 
#     n_repeats=1000, 
#     domains=[{'start':300, 'end':350}], 
# )
