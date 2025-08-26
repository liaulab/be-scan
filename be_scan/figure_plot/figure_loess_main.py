
from be_scan.figure_plot.figure_classes import *
from be_scan.figure_plot.figure_loess_helpers import *

def loess_smoothing(
    df, x_column, comparisons, span,
    interp_method='quadratic', n_repeats=1000, random_seed=None,
    loess_kws=None, interp_kws=None, smm_multipletests_kws=None,
):
    # OR OPERATOR ASSIGNS VALUE IF NOT NONE, ELSE ASSIGNS A DEFAULT #
    loess_kws = loess_kws or {
        'missing':'raise', 'return_sorted':False, 'it':0
        }
    interp_kws = interp_kws or {
        'fill_value':'extrapolate'
        }
    smm_multipletests_kws = smm_multipletests_kws or {
        'alpha': 0.05, 'method': 'fdr_bh', 'is_sorted': False, 'returnsorted': False
        }

    # CHECKS #
    assert x_column in df.columns, f'{x_column} not in dataframe'
    for comp in comparisons:
        assert comp in df.columns, f'{comp} not in dataframe'
    assert span > 0, f'span must be greater than 0'
    assert n_repeats > 0, f'n_repeats must be greater than 0'
    assert interp_method in interp_methods, f'interp_method must be one of {interp_methods}'

    # FILTER AND FLOOR DATA #
    df_filtered = df[[x_column] + comparisons].copy()
    df_filtered = df_filtered[df_filtered[x_column] >= 0].sort_values(by=x_column)

    results = {}
    for comp in comparisons:
        x_obs, y_obs = df_filtered[x_column], df_filtered[comp]
        x_min, x_max = np.floor(x_obs.min()), np.ceil(x_obs.max())
        span_new = span / (x_max - x_min) if span > 1 else span
        print(f"Span (decimal):", span)

        df_loess = loess_v3(x_obs, y_obs, span=span_new,
                            interp_method=interp_method, loess_kws=loess_kws, interp_kws=interp_kws)
        df_rand = randomize(x_obs, y_obs, span=span_new, n_repeats=n_repeats,
                            random_seed=random_seed, interp_method=interp_method)
        df_pvals = calculate_sig(df_loess, df_rand, n_repeats=n_repeats,
                                 smm_multipletests_kws=smm_multipletests_kws)

        results[comp] = {
            'loess': df_loess, 'rand': df_rand, 'pvals': df_pvals
        }
    return results
