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

import warnings
warnings.filterwarnings("ignore")

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

# specific to the loess function
import scipy.interpolate as interp
from statsmodels.nonparametric.smoothers_lowess import lowess
import statsmodels.stats.multitest as smm
import random

def loess_smoothing(df_filepath, 
                    x_column, comparisons, span, 

    x_out=None, interp_method='quadratic', 
    loess_kws={'missing':'raise', 'return_sorted':False, 'it':0}, 
    interp_kws={'fill_value':'extrapolate'}, 
    n_repeats=10000, 
    smm_multipletests_kws={'alpha':0.05, 'method':'fdr_bh', 
                           'is_sorted':False, 'returnsorted':False}, 

    savefig=True, out_name='scatterplot', out_type='png', out_directory='', show=True, # output params
    return_df=True, 
    ): 

    """[Summary]
    this function calculates the smoothed arbitrary values for each value 
    in your x_out for the enrichment profile of your POI
    
    Parameters
    ----------
    df_filepath : str, required
        filepath to .csv data
    x_column : str, required
        column of .csv, typically amino acid position
    comparisons : list of str, required
        list of comparisons that correspond to columns of data
    span: float, required (recommended range to start is ~0.05-0.10)
        hyperparameter for the span that you want to smooth over, requires optimization

    x_out: list or numpy array, optional, defaults to None
        the values you want to interpolate over to get values for, 
        usually an array of discrete values for the whole protein (i.e. np.arange(1, len(POI)+1))
    interp_method: str, optional, defaults to 'quadratic'
        method should 'statsmodels' or on be one of 
        'linear', 'nearest', 'nearest-up', 'zero', 
        'slinear', 'quadratic', 'cubic', 'previous', or 'next'
        from https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html

    loess_kws: dict, optional, defaults to
        {'missing':'raise', 'return_sorted':False, 'it':0}
        input params for statsmodels.nonparametric.smoothers_lowess.lowess()
        https://www.statsmodels.org/devel/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html
    interp_kws: dict, optional, defaults to
        {'fill_value':'extrapolate'}
        input params for scipy.interpolate.interp1d()
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html

    n_repeats: int, optional, defaults to 10,000
        the number of times you want to shuffle
        10000x is recommended for final data but for preliminaty sample 4000x is the recommended minimum
    smm_multipletests_kws: dict, optional, defaults to 
        {'alpha':0.05, 'method':'fdr_bh', 'is_sorted':False, 'returnsorted':False}
        input params for smm.multipletests()
        https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
    
    savefig: bool, optional, defaults to True
        whether or not to save the figure
    out_name : str, optional, defaults to 'loess_smoothed'
        name of figure output
    out_type : str, optional, defaults to 'pdf'
        file type of figure output
    out_directory : str, optional, defaults to ''
        path to output directory
    show : bool, optional, defaults to True
        whether or not to show the plot
    return_df: bool, optional, defaults to True
        whether or not to return the dataframes

    Returns
    ----------

    """

    mpl.rcParams.update({'font.size': 10})
    fig, ax = plt.subplots(nrows=len(comparisons), ncols=1, 
                          figsize=(10, 5*len(comparisons)))
    
    # check variable inputs
    assert interp_method in ['linear', 'nearest', 'nearest-up', 'zero', 
                             'slinear', 'quadratic', 'cubic', 'previous', 'next']

    # process input dataframe and keep only necessary data fom the input
    df_filepath = Path(df_filepath)
    df_data = pd.read_csv(df_filepath)

    # check all column names are valid, and filter down data
    assert x_column in df_data.columns, "check param: xlim_kws"
    for comp in comparisons: 
        assert comp in df_data.columns, f"check param: comparisons {comp}"
    df_filtered = df_data[[x_column] + comparisons]
    df_filtered = df_filtered[df_filtered[x_column] >= 0].sort_values(by=x_column)

    return_list = []
    # process columns from dataframe and input into loess_v3
    if len(comparisons) == 1: 
        ax = [ax]
    for i, comp in enumerate(comparisons): 
        x_obs, y_obs = df_filtered[x_column], df_filtered[comp]

        df_loess = loess_v3(x_obs=x_obs, y_obs=y_obs, span=span, 
                            x_out=x_out, interp_method=interp_method, 
                            loess_kws=loess_kws, interp_kws=interp_kws, 
                            ) 
        df_loess_rand = randomize(x_obs=x_obs, y_obs=y_obs, span=span, 
                                  x_out=x_out, n_repeats=n_repeats, interp_method=interp_method, 
                                  )
        df_pvals = calculate_sig(df_loess=df_loess, df_loess_rand=df_loess_rand, 
                                 n_repeats=n_repeats, smm_multipletests_kws=smm_multipletests_kws, 
                                 )

        df_plotting = df_pvals[['corr_pval']].copy()
        # plot -log10 of p-values so that hits will be on top
        # the last part is added because some of my p-values were 0 and cannot be log-transformed
        df_plotting['-log10'] = np.log10(df_plotting['corr_pval'] + 10**-4) * -1
        ax[i].plot(df_plotting.index, df_plotting['-log10'], color='gainsboro', ) 
        ax[i].axhline(y=np.log10(0.05 + 10**-4)*-1,ls='--', c='k', linewidth=1) 

        fig.tight_layout()
        ### title and labels

        if return_df: 
            return_list.append((df_loess, df_loess_rand, df_pvals))

    path = Path.cwd()
    if savefig:
        outpath = path / out_directory
        out = out_name + comp + '.' + out_type
        plt.savefig(outpath / out, format=out_type)
    if show: 
        plt.show()
    plt.close()
    
    if return_df: 
        return return_list
    else: 
        return None

def loess_v3(x_obs, y_obs, span, 
             
    x_out=None, interp_method='quadratic', 
    loess_kws={}, interp_kws={}, 
    ):

    """[Summary]
    this function calculates the smoothed arbitrary values for each value 
    in your x_out for the enrichment profile of your POI
    
    Parameters
    ----------
    x_obs: list or numpy array, required
        the "x-axis" positions of your data, usually some summarized aa number of edited window
    y_obs: list or numpy array, required
        the "y-axis" values of your data, usually normalized lfc of gRNA enrichment
    span: float, required (recommended range to start is ~0.05-0.10)
        hyperparameter for the span that you want to smooth over, requires optimization

    x_out: list or numpy array, optional, defaults to None
        the values you want to interpolate over to get values for, 
        usually an array of discrete values for the whole protein (i.e. np.arange(1, len(POI)+1))
    interp_method: str, optional, defaults to 'quadratic'
        method should 'statsmodels' or on be one of 
        'linear', 'nearest', 'nearest-up', 'zero', 
        'slinear', 'quadratic', 'cubic', 'previous', or 'next'
        from https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html

    loess_kws: dict, optional, defaults to
        {'missing':'raise', 'return_sorted':False, 'it':0}
        input params for statsmodels.nonparametric.smoothers_lowess.lowess()
        https://www.statsmodels.org/devel/generated/statsmodels.nonparametric.smoothers_lowess.lowess.html
    interp_kws: dict, optional, defaults to
        {'fill_value':'extrapolate'}
        input params for scipy.interpolate.interp1d()
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html

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
    # df_loess = df_loess.reset_index(drop=True)

    if interp_method == 'statsmodels':
        df_loess['y_loess'] = lowess(endog=y_obs, exog=x_obs, 
                                     frac=span, **loess_kws, 
                                     )
    else: 
        loess_out = lowess(endog=y_obs, exog=x_obs, 
                           frac=span, **loess_kws)
            # for your x values, do loess smoothing and get y values

        # for some reason sometimes loess_out is 2D and sometimes its 1D
        loess_out = np.array(loess_out)
        df_interp = pd.DataFrame(loess_out)
        if df_interp.shape[1] == 2: 
            df_interp = df_interp.drop(columns=[0])
            df_interp = df_interp.rename({1:0}, axis='columns')
        df_interp['x_obs'] = x_obs.values
        df_interp = df_interp.rename({0: 'y_obs_loess'}, axis='columns')
        df_interp = df_interp.groupby('x_obs', as_index=False).agg({'y_obs_loess':'mean'}) 
            # in case theres duplicate values in x_obs
        fx_interp = interp.interp1d(x=df_interp['x_obs'], 
                                    y=df_interp['y_obs_loess'], 
                                    kind=interp_method, **interp_kws)
            # interpolate missing variables
        df_loess['y_loess'] = fx_interp(df_loess['x_vals'])
            # apply loess interpolated function to a range of x values

    df_loess['type'] = np.where(df_loess['x_vals'].isin(x_obs), 'loess', 'interp')
    return df_loess

def randomize(x_obs, y_obs, span, 
    
    x_out=None, n_repeats=10000, interp_method='quadratic', 
    ): 

    """[Summary]
    this function provides the null distribution to compare your actual loess again; 
    it calculates the smoothed arbitrary values for each value 
    in your x_out (see below) for the shuffled enrichment profile of your POI
    
    Parameters
    ----------
    x_obs: str, required
        the "x-axis" position of your data, usually some summarized aa number of edited window
    y_obs: str, required
        the "y-axis" value of your data, usually normalized lfc of gRNA enrichment
    span: float, required (recommended range to start is ~0.05-0.10)
        hyperparameter for the span that you want to smooth over, requires optimization

    x_out: list or numpy array, optional, defaults to None
        the values you want to interpolate over to get values for, 
        usually an array of discrete values for the whole protein (i.e. np.arange(1, len(POI)+1))
    n_repeats: int, optional, defaults to 10,000
        the number of times you want to shuffle
        10000x is recommended for final data but for preliminaty sample 4000x is the recommended minimum
    interp_method: str, optional, defaults to 'quadratic'
        method should 'statsmodels' or on be one of 
        'linear', 'nearest', 'nearest-up', 'zero', 
        'slinear', 'quadratic', 'cubic', 'previous', or 'next'
        from https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
        
    Returns
    ----------
    df_loess_rand: pandas dataframe
        has n_repeats number of columns each of calculated loess scores from randomized y_obs
    """

    df_loess_rand = pd.DataFrame()
    df_loess_rand = pd.concat([loess_v3(x_obs=x_obs, y_obs=random.sample(list(y_obs), len(y_obs)), 
                                        span=span, interp_method=interp_method).x_vals] + 
                              [loess_v3(x_obs=x_obs, y_obs=random.sample(list(y_obs), len(y_obs)), 
                                        span=span, interp_method=interp_method).y_loess for _ in range(0, n_repeats)] +
                              [loess_v3(x_obs=x_obs, y_obs=random.sample(list(y_obs), len(y_obs)), 
                                        span=span, interp_method=interp_method).type], 
                               ignore_index=True, axis=1)
    return df_loess_rand

def calculate_sig(df_loess, df_loess_rand, 

    n_repeats=10000, smm_multipletests_kws={}
    ):

    """[Summary]
    this function calculates the significance of a "spike" in your loess-ed enrichment profile
    by calculating if it's > 95% of null distribution
    
    Parameters
    ----------
    df_loess: pandas df, required
        output from loess_v3()
    df_loess_rand: pandas df, required
        output from randomize()

    n_repeats: int, optional, defaults to 10,000
        the number of times you want to shuffle
        10000x is recommended for final data but for preliminaty sample 4000x is the recommended minimum
    smm_multipletests_kws: dict, optional, defaults to 
        {'alpha':0.05, 'method':'fdr_bh', 'is_sorted':False, 'returnsorted':False}
        input params for smm.multipletests()
        https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html

    Returns
    ----------
    df_pvals: pandas dataframe
        dataframe with your position values (x_out) and corresponding loess values (y_loess) from input, 
        and empirical p-val (sig) as well as the corrected p-values (corr_pval)
    """

    df_pvals = pd.DataFrame({'x_vals': list(df_loess.x_vals), 'y_loess': list(df_loess.y_loess)}, 
                            columns=['x_vals', 'y_loess'])
    df_pvals['obs_gt'] = df_loess_rand[df_loess_rand.columns[1:-1]].gt(df_pvals['y_loess'], 
                                                                       axis=0).sum(axis=1) 
        # get # of values greater than y_loess
    df_pvals['1t_pval'] = df_pvals['obs_gt'] / n_repeats 
        # divide "rank" of obs val by N to get empirical p-val

    temp = smm.multipletests(df_pvals['1t_pval'], **smm_multipletests_kws) 
        # apply benjamini-hochberg FDR correction
    df_pvals['sig'] = temp[0]
    df_pvals['corr_pval'] = temp[1]
    return df_pvals
