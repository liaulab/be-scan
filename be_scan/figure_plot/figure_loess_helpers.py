
from be_scan.figure_plot.figure_classes import *

def loess_v3(
    x_obs, y_obs, span,
    interp_method,
    loess_kws={}, interp_kws={},
    ):

    """[Summary]
    This function calculates the smoothed arbitrary values for each value
    in your x_out for the enrichment profile of your POI

    Returns
    ----------
    df_loess: pandas dataframe
        dataframe of x_vals and y_obs from input,
        and y_loess and type from output corresponding to
        the loess calculated or interpolated values of the loess smoothing and the type
    """

    x_obs = x_obs.astype(float)
    # RESOLUTION OF INTEGERS ONLY #
    df_loess = pd.DataFrame()
    df_loess['x_vals'] = np.arange(math.floor(x_obs.min()),
                                   math.ceil(x_obs.max())+1, 1)

    # RUN LOWESS #
    loess_out = lowess(endog=y_obs, exog=x_obs, frac=span, **loess_kws)
    if interp_method == 'statsmodels':
        df_loess['y_loess'] = loess_out

    else:
        # FIXING 1D vs 2D ARRAY #
        loess_out = np.array(loess_out)
        df_interp = pd.DataFrame(loess_out)
        if df_interp.shape[1] == 2:
            df_interp = df_interp.drop(columns=[0])
            df_interp = df_interp.rename({1:0}, axis='columns')
        # X AND Y VALUES #
        df_interp['x_obs'] = x_obs.values
        df_interp = df_interp.rename({0: 'y_obs_loess'}, axis='columns')
        # IN CASE OF DUP VALUES FOR X #
        df_interp = df_interp.groupby('x_obs', as_index=False).agg({'y_obs_loess':'mean'})
        # INTERPOLATE MISSING VARIABLES #
        fx_interp = interp.interp1d(x=df_interp['x_obs'], y=df_interp['y_obs_loess'],
                                    kind=interp_method, **interp_kws)
        # APPLY LOESS INTERPOLATED FUNCTION #
        df_loess['y_loess'] = fx_interp(df_loess['x_vals'])

    # LOESS OR INTERPOLATED #
    df_loess['type'] = np.where(df_loess['x_vals'].isin(x_obs), 'loess', 'interp')
    return df_loess

def randomize(
    x_obs, y_obs, span,
    n_repeats=1000, interp_method='quadratic',
    random_seed=None,
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

    random.seed(random_seed)
    loess_args = {'x_obs':x_obs, 'span':span, 'interp_method':interp_method}

    loess_rand = []
    for _ in range(n_repeats):
        loess_args['y_obs'] = random.sample(list(y_obs), len(y_obs))
        loess_out = loess_v3(**loess_args)
        loess_rand.append(loess_out.y_loess)

    df_rand = pd.concat([loess_out.x_vals] + loess_rand + [loess_out.type],
                        ignore_index=True, axis=1 )
    return df_rand

def calculate_sig(
    df_loess, df_rand,
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
    # NUMBER OF VALUES GREATER THAN y_loess #
    gt_counts = np.sum(df_rand[df_rand.columns[1:-1]] > df_loess['y_loess'].values[:, None], axis=1)
    # DIviDE "RANK" OF OBSERVED VALUE BY N TO GET EMPIRICAL P-VAL #
    df_pvals['1t_pval'] = gt_counts / n_repeats

    # APPLY BENJAMINI-HOCHBERG FDR CORRECTION #
    temp = smm.multipletests(df_pvals['1t_pval'], **smm_multipletests_kws)
    df_pvals['sig'] = temp[0]
    df_pvals['corr_pval'] = temp[1]
    return df_pvals

def increase_resolution(
    x, y, thr, res=0.01
    ):

    """[Summary]
    Helper function to increase the resolution of a line above and
    below a threshold. This helps with coloring.
    """

    new_x, new_y = [], []
    for i in range(len(x) - 1):
        new_x.append(x[i])
        new_y.append(y[i])

        # y=mx+b OR y=-mx+b #
        m = (y[i+1] - y[i]) / (x[i+1] - x[i])
        b = y[i] - m * x[i]
        # ADD POINTS TO OUTPUT #
        if y[i] < thr < y[i+1]:
            new_x.extend([(thr-res-b)/m, (thr-b)/m, (thr+res-b)/m])
            new_y.extend([thr-res, thr, thr+res])
        if y[i] > thr > y[i+1]:
            new_x.extend([(thr+res-b)/m, (thr-b)/m, (thr-res-b)/m])
            new_y.extend([thr+res, thr, thr-res])

    new_x.append(x[-1])
    new_y.append(y[-1])
    return new_x, new_y

def split_on_value(
    lst_x, lst_y, val
    ):

    """[Summary]
    Helper function to split line segments that color the regression line
    crossing a threshold differently. This helps with coloring.
    """

    result_x, result_y = [], []
    current_x, current_y = [], []

    for item_x, item_y in zip(lst_x, lst_y):
        current_x.append(item_x)
        current_y.append(item_y)
        if item_y == val:
            result_x.append(current_x[:-1])
            result_y.append(current_y[:-1])
            current_x = []
            current_y = []
    if len(current_y) > 1 or (current_y and current_y != [val]):
        result_x.append(current_x)
        result_y.append(current_y)
    return result_x, result_y

interp_methods = ['linear', 'nearest', 'nearest-up', 'zero',
                  'slinear', 'quadratic', 'cubic', 'previous', 'next']
