"""
Author: Calvin XiaoYang Hu, Simon Shen, Kevin Ngan
Adapted from: Kevin Ngan from KCN_masterfunctions_v6_200406.py
Date: 231128

{Description: Aggregate raw counts and perform log2-transform and t0 normalization}
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma
import statsmodels.stats.multitest as smm

def calc_pvals(df_reps):
    """
    Two-sided one-sample t-test of each row's replicate values against zero.

    Parameters
    ------------
    df_reps : pandas DataFrame
        One column per replicate, one row per sgRNA.

    Returns
    ------------
    numpy array of p-values, one per row. Rows with fewer than two usable
    replicates, or with no variance across them, yield NaN rather than a
    misleading value.

    Notes
    ------------
    Zero-variance rows are deliberately NaN rather than 0. scipy reports p=0
    for identical replicates, because the t-statistic is infinite, and -log10
    of that is ~308 -- the most significant value possible. An sgRNA with zero
    counts in every replicate has exactly identical log-fold-changes, so
    dropped-out guides would otherwise dominate every volcano plot. With no
    spread between replicates there is no evidence about significance, and NaN
    is the honest answer.
    """
    values = df_reps.to_numpy(dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        _, pvals = stats.ttest_1samp(values, popmean=0.0, axis=1,
                                     nan_policy='omit')
    pvals = np.asarray(pvals, dtype=float)

    # mask rows that carry no information: too few replicates, or no spread
    n_obs = np.sum(~np.isnan(values), axis=1)
    with np.errstate(invalid='ignore'):
        spread = np.nanstd(values, axis=1)
    pvals = np.where((n_obs < 2) | (spread == 0), np.nan, pvals)
    return pvals


def _trigamma_inverse(y):
    """
    Solve trigamma(x) = y by Newton iteration.

    Used when fitting the prior variance distribution. Mirrors limma's
    trigammaInverse; the initial guess and update are Smyth's.
    """
    y = float(y)
    if not np.isfinite(y) or y <= 0:
        return np.inf
    x = 0.5 + 1.0 / y
    for _ in range(50):
        tri = polygamma(1, x)
        step = tri * (1 - tri / y) / polygamma(2, x)
        x = x + step
        if abs(step / x) < 1e-8:
            break
    return x


def fit_variance_prior(variances, df):
    """
    Fit a scaled inverse chi-squared prior to a set of sample variances.

    This is the method-of-moments fit from Smyth (2004), the step that makes a
    moderated t-test possible: it learns how variable sgRNA variances are
    across the whole screen, and uses that as a prior.

    Parameters
    ------------
    variances : array of per-sgRNA sample variances
    df : int, residual degrees of freedom per sgRNA (n_replicates - 1)

    Returns
    ------------
    (prior_var, prior_df). prior_df is inf when the observed variances are no
    more dispersed than sampling noise alone would explain, which collapses
    the moderated test to using a single pooled variance for every guide.
    """
    variances = np.asarray(variances, dtype=float)
    usable = variances[np.isfinite(variances) & (variances > 0)]
    if usable.size < 2 or df <= 0:
        return np.nan, np.inf

    z = np.log(usable)
    e = z - digamma(df / 2.0) + np.log(df / 2.0)
    e_mean = float(np.mean(e))
    e_var = float(np.var(e, ddof=1)) - float(polygamma(1, df / 2.0))

    if e_var > 0:
        prior_df = 2.0 * _trigamma_inverse(e_var)
        prior_var = float(np.exp(e_mean + digamma(prior_df / 2.0)
                                 - np.log(prior_df / 2.0)))
    else:
        prior_df = np.inf
        prior_var = float(np.exp(e_mean))
    return prior_var, prior_df


def calc_moderated_pvals(df_reps):
    """
    Moderated (variance-shrunk) two-sided t-test against zero.

    With three replicates the ordinary t-test is dominated by how lucky each
    sgRNA's variance estimate was: a guide whose replicates happen to agree
    closely outranks one with a far larger effect. Screens routinely surface
    non-targeting controls this way.

    The moderated test borrows strength across the whole library. Each guide's
    variance is pulled toward a prior fitted from every other guide, so an
    implausibly small variance cannot by itself produce a tiny p-value, and
    the degrees of freedom rise accordingly.

    Unlike the ordinary test, this is defined for a guide whose replicates are
    identical: shrinkage gives it the prior variance rather than zero.

    Returns
    ------------
    numpy array of p-values, one per row, NaN where fewer than two replicates
    are usable.
    """
    values = df_reps.to_numpy(dtype=float)
    n_obs = np.sum(~np.isnan(values), axis=1)

    with np.errstate(invalid='ignore', divide='ignore'):
        means = np.nanmean(values, axis=1)
        variances = np.nanvar(values, ddof=1, axis=1)

    # fit the prior on the guides that carry a usable variance estimate
    full = int(np.max(n_obs)) if n_obs.size else 0
    df_resid = full - 1
    prior_var, prior_df = fit_variance_prior(variances[n_obs == full], df_resid)

    if not np.isfinite(prior_var):
        return np.full(values.shape[0], np.nan)

    df_obs = np.maximum(n_obs - 1, 0).astype(float)
    if np.isinf(prior_df):
        posterior_var = np.full_like(variances, prior_var)
        total_df = np.full_like(df_obs, np.inf)
    else:
        posterior_var = ((prior_df * prior_var) + (df_obs * np.nan_to_num(variances))) \
            / (prior_df + df_obs)
        total_df = prior_df + df_obs

    with np.errstate(invalid='ignore', divide='ignore'):
        t_mod = means / np.sqrt(posterior_var / n_obs)
        pvals = 2.0 * stats.t.sf(np.abs(t_mod), total_df)

    return np.where(n_obs >= 2, pvals, np.nan)


def fdr_correct(pvals):
    """
    Benjamini-Hochberg false discovery rate correction.

    A screen tests every sgRNA at once, so raw p-values below 0.05 are
    expected by chance in proportion to the library size. NaN entries are
    excluded from the correction and stay NaN.
    """
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(pvals.shape, np.nan)
    valid = np.isfinite(pvals)
    if valid.sum() == 0:
        return out
    _, qvals, _, _ = smm.multipletests(pvals[valid], method='fdr_bh')
    out[valid] = qvals
    return out


def neg_log10(pvals):
    """
    -log10 of a p-value array, for volcano-style plotting.

    A p-value of exactly 0 would give infinity, so it is floored at the
    smallest positive double first. NaN stays NaN.
    """
    pvals = np.asarray(pvals, dtype=float)
    floored = np.where(pvals > 0, pvals, np.finfo(float).tiny)
    with np.errstate(invalid='ignore', divide='ignore'):
        out = -np.log10(floored)
    return np.where(np.isnan(pvals), np.nan, out)


def log_transform(
    sample_sheet, 
    library_counts, 
                   
    controls=['t0'], 
    in_dir='', out_dir='', out_file='library_LFC.csv', 
    save=True, return_df=True,
    ):
    
    """[Summary]
    For a given set of samples and their raw read count files from count_reads,
    aggregate them into a single dataframe, normalize counts to reads per million,
    perform log2 transform, and then normalize to t0. The aggregated raw reads,
    log2 transformed values, and t0 normalized values can be saved to csv.

    Parameters
    ------------
    sample_sheet : str or path or pandas DataFrame()
        REQUIRED COLS: 'condition', 'counts_file'
        a sheet with information on sequence id, 
        in_fastq (string or path to the FASTQ file to be processed), 
        out_counts (string or path for the output csv file with perfect sgRNA matches ex: 'counts.csv'),
        out_np (string or path for the output csv file with non-perfect sgRNA matches ex: 'noncounts.csv'), 
        out_stats (string or path for the output txt file with the read counting statistics ex: 'stats.txt'), 
        condition names, and condition categories
    library_counts : str or path or pandas DataFrame()
        String or path to the reference file. library_counts must have column headers,
        with 'sgRNA_seq' as the header for the column with the sgRNA sequences.

    controls : str, default ['t0']
        Name of the control condition samples in sample_sheet. 
    in_dir : str or path, defaults to ''
        String or path to the directory where library_counts and sample_sheet is found. 
    out_dir : str or path, defaults to ''
        String or path to the directory where all files are found. 
    out_file : str or path, defaults to 'agg_log2_t0.csv'
        Name of output dataframe with guides and counts. 
    return_df : bool, default True
        Whether or not to return the resulting dataframe
    save : bool, default True
        Whether or not to save the resulting dataframe

    Returns
    ------------
    df_ref : pandas DataFrame

        For every condition with more than one replicate, six columns are added
        alongside the existing '<cond>_avg' and '<cond>_stdev'::

            <cond>_pval          two-sided one-sample t-test of that sgRNA's
                                 replicate values against 0
            <cond>_neglogp       -log10 of the above
            <cond>_pval_mod      the same test with the variance shrunk toward a
                                 prior fitted across the whole library
            <cond>_neglogp_mod   -log10 of the above
            <cond>_qval          Benjamini-Hochberg FDR of <cond>_pval
            <cond>_qval_mod      Benjamini-Hochberg FDR of <cond>_pval_mod

        Which to use: with only two or three replicates, prefer the moderated
        columns. The ordinary t-test ranks largely on how lucky each guide's
        variance estimate was, so a guide whose replicates happen to agree closely
        can outrank one with a far larger effect, and non-targeting controls
        routinely appear near the top. Read either alongside '<cond>_avg' -- a
        small p-value on a negligible effect is still a negligible effect.

        Conditions with a single replicate get none of these, since a p-value
        across one measurement is not defined. No existing column is renamed,
        reordered or altered.
    """
    in_path = Path(in_dir)
    out_path = Path(out_dir)

    # IMPORT FILES, CHECK FOR REQUIREMENTS #
    if isinstance(library_counts, (str, Path)):
        df_ref = pd.read_csv(in_path / library_counts)
    else: df_ref = library_counts.copy()
    if isinstance(sample_sheet, (str, Path)):
        df_samples = pd.read_csv(in_path / sample_sheet)
    else: df_samples = sample_sheet.copy()

    dict_counts = dict(zip(df_samples.condition, df_samples.counts_file))
    dict_conds = dict(zip(df_samples.condition, df_samples.agg_conditions))
    df_map = pd.DataFrame(data=dict_conds.items(), columns=['rep','condition'])

    # CHECK CONTROLS IN DF AND AVERAGE IF NECESSARY #
    conditions = df_samples.condition.tolist()
    for control in controls: 
        if control not in dict_counts.keys(): 
            raise Exception (f"{sample_sheet} is missing the {control} sample")
    if len(controls) > 0: 
        df_ref['control_avg'] = df_ref[controls].mean(axis=1)
        conditions = ['control_avg'] + conditions
    
    # CALCULATE LFC AND CONTROLS #
    for sample in conditions: 
        # log2 normalization
        total_reads = pd.to_numeric(df_ref[sample]).sum()
        suffix = '_LFC'
        df_ref[sample+suffix] = df_ref[sample].apply(lambda x: np.log2((x * 1000000 / total_reads) + 1))
        if sample not in controls and len(controls) > 0: 
            # t0 normalization
            suffix = '_LFCminusControl'
            df_ref[sample+suffix] = df_ref[sample+'_LFC'].sub(df_ref['control_avg_LFC'])

    # AVERAGE LFC #
    for cond in df_map['condition'].unique().tolist(): 
        if cond not in controls: 
            # for each unique condition, find the reps
            reps = [x+suffix for x in df_map.loc[df_map['condition'] == cond]['rep'].tolist()]
            # skip averaging for single replicates (otherwise breaks script)
            if len(reps) > 1:
                df_ref[cond+suffix+'_avg'] = df_ref[reps].mean(axis=1)
                df_ref[cond+suffix+'_stdev'] = df_ref[reps].std(axis=1)
                # SIGNIFICANCE ACROSS REPLICATES #
                # Two-sided one-sample t-test of each sgRNA's replicate values
                # against 0, i.e. "is this guide's effect distinguishable from
                # no change?". Only defined with more than one replicate, so
                # single-rep conditions get no p-value columns at all.
                pvals = calc_pvals(df_ref[reps])
                df_ref[cond+suffix+'_pval'] = pvals
                df_ref[cond+suffix+'_neglogp'] = neg_log10(pvals)
                # MODERATED TEST AND MULTIPLE-TESTING CORRECTION #
                # With few replicates the ordinary t-test ranks on how lucky a
                # guide's variance estimate was, so a variance-shrunk version
                # is provided alongside it. Both are then FDR-corrected across
                # the library, since every sgRNA is tested at once.
                pvals_mod = calc_moderated_pvals(df_ref[reps])
                df_ref[cond+suffix+'_pval_mod'] = pvals_mod
                df_ref[cond+suffix+'_neglogp_mod'] = neg_log10(pvals_mod)
                df_ref[cond+suffix+'_qval'] = fdr_correct(pvals)
                df_ref[cond+suffix+'_qval_mod'] = fdr_correct(pvals_mod)
            elif len(reps) == 1:
                df_ref[cond+suffix+'_avg'] = df_ref[reps]
                df_ref[cond+suffix+'_stdev'] = 0
            else:
                raise Exception('Error! Replicate number not valid')

    # SAVE DF AND RETURN #
    if save == True: 
        Path.mkdir(out_path, exist_ok=True)
        df_ref.to_csv(out_path / out_file, index=False)
        print('merge_and_norm output to', str(out_path / out_file))
    print('Merge and normalize completed')
    if return_df: 
        return df_ref
    