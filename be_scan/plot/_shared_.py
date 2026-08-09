"""
Author: Calvin XiaoYang Hu

{Description: shared building blocks used by every be_scan.plot function}

Before 3.0.0 each plotting function implemented whichever of these it happened
to need, so scatterplot could filter and colour and annotate while corr_heatmap
could do none of those, and two functions that took the same parameter did not
always give it the same default. Routing every function through the helpers
here is what makes the plots behave consistently.

Every helper is a no-op when its feature is switched off, so adding a call to
one of these to an existing function cannot change that function's output
unless the caller opts in.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import warnings

from be_scan.plot._annotating_ import calc_neg_ctrls, norm_to_intergenic_ctrls


def load_dataframe(df_filepath):
    """
    Accept either a DataFrame or a path to a .csv, and return a copy.

    A copy, so no plotting function ever mutates the caller's frame.
    """
    if isinstance(df_filepath, pd.DataFrame):
        return df_filepath.copy()
    if isinstance(df_filepath, (str, Path)):
        return pd.read_csv(Path(df_filepath))
    raise ValueError("df_filepath must be a pandas DataFrame or a path string/Path.")


def normalize_to_controls(df_data, comparisons, neg_ctrl, neg_ctrl_col,
                          neg_ctrl_conditions):
    """
    Centre each comparison on its negative controls.

    Returns (df_data, list_ctrlstats). list_ctrlstats is None when neg_ctrl is
    False, which is what callers check before drawing +/- 2 SD guide lines.
    """
    if not neg_ctrl:
        return df_data, None

    assert isinstance(neg_ctrl_col, str) and neg_ctrl_col in df_data.columns.tolist(), \
        "check param: neg_ctrl_col"
    assert isinstance(neg_ctrl_conditions, list), "check param: neg_ctrl_conditions"

    _, list_ctrlstats, avg_dict = calc_neg_ctrls(
        df_data, comparisons, neg_ctrl_col, neg_ctrl_conditions)
    df_data = norm_to_intergenic_ctrls(df_data, comparisons, avg_dict)
    return df_data, list_ctrlstats


def control_stdev(list_ctrlstats, comparison):
    """The negative-control standard deviation for one comparison, or None."""
    if not list_ctrlstats:
        return None
    for tup in list_ctrlstats:
        if tup[0] == comparison:
            return tup[2]
    return None


def apply_xwindow(df_data, x_column, xwindow):
    """Restrict rows to an inclusive range on x_column."""
    if xwindow is None:
        return df_data
    if len(xwindow) != 2:
        warnings.warn("Warning: xwindow should be a list of 2 values", UserWarning)
        return df_data
    return df_data[(df_data[x_column] >= xwindow[0]) & (df_data[x_column] <= xwindow[1])]


def apply_filters(df_data, filter_val=False, val_cols=(), val_min=0.0,
                  filter_params=False, params_cols=(), params_conditions=(),
                  verbose=True):
    """
    Drop rows by a minimum value and/or by categorical membership.

    filter_val keeps rows where every column in val_cols exceeds val_min.
    filter_params keeps rows whose value in each params_cols entry is in the
    matching params_conditions list.
    """
    if filter_val:
        assert isinstance(val_min, (int, float)), "check param: val_min"
        assert isinstance(val_cols, (list, tuple)) and len(val_cols) > 0, \
            "check param: val_cols"
        for col in val_cols:
            if verbose:
                print(f'Filter data by {col} to be greater than {str(val_min)}')
            df_data = df_data[df_data[col] > val_min]

    if filter_params:
        assert isinstance(params_cols, (list, tuple)), "check param: params_cols"
        assert isinstance(params_conditions, (list, tuple)), "check param: params_conditions"
        for col, conds in zip(params_cols, params_conditions):
            if verbose:
                print(f'Filter data by {col} to contain:')
                print(conds)
            df_data = df_data.loc[df_data[col].isin(conds)]

    return df_data


def resolve_hue(df_data, include_hue, hue_col, pal):
    """
    Colour settings for a seaborn call.

    Returns a dict to splat into the plotting call: empty when include_hue is
    False, otherwise hue, hue_order and palette.
    """
    if not include_hue:
        return {}
    assert isinstance(hue_col, str) and hue_col in df_data.columns, "check param: hue_col"
    unique_types = df_data[hue_col].dropna().unique()
    return {'hue': hue_col, 'hue_order': unique_types,
            'palette': sns.color_palette(pal, len(unique_types))}


def hue_color_map(df_data, include_hue, hue_col, pal):
    """A value -> hex colour mapping, for backends that need explicit colours."""
    import matplotlib.colors as mcolors
    if not include_hue:
        return None
    assert isinstance(hue_col, str) and hue_col in df_data.columns.tolist(), \
        "check param: hue_col"
    unique_types = df_data[hue_col].dropna().unique()
    colors = sns.color_palette(pal, n_colors=len(unique_types))
    return dict(zip(unique_types, [mcolors.to_hex(c) for c in colors]))


def split_highlight(df_data, highlight, highlight_col, highlight_vals):
    """
    Split rows into (highlighted, rest).

    Used to draw a subset of guides with heavier styling on top of the others.
    Returns (None, df_data) when highlighting is off.
    """
    if not highlight:
        return None, df_data
    assert isinstance(highlight_col, str) and highlight_col in df_data.columns, \
        "check param: highlight_col"
    assert isinstance(highlight_vals, (list, tuple)), "check param: highlight_vals"
    mask = df_data[highlight_col].isin(highlight_vals)
    return df_data[mask], df_data[~mask]


def draw_domains(ax, domains, domains_alpha=0.25, domains_color='lightblue'):
    """Shade protein domains as vertical spans."""
    for d in domains or []:
        ax.axvspan(d['start'], d['end'],
                   alpha=d.get('alpha', domains_alpha),
                   facecolor=d.get('color', domains_color))


def draw_control_lines(ax, stdev, n_sd=2, **axhline_kws):
    """Horizontal guide lines at +/- n_sd times the control standard deviation."""
    if stdev is None:
        return
    ax.axhline(y=n_sd * stdev, **axhline_kws)
    ax.axhline(y=-n_sd * stdev, **axhline_kws)


def select_annotated(df_data, comparison, annot_top=None, annot_cutoff=None,
                     annot_abs=None):
    """
    The rows a plot should label.

    annot_cutoff labels everything whose absolute value exceeds it; annot_abs
    labels the n largest by absolute value; annot_top labels the n highest, or
    the n lowest when negative. Returns an empty frame when none is set.
    """
    if annot_cutoff is not None:
        return df_data[df_data[comparison].abs() > annot_cutoff]
    if annot_abs is not None:
        ordered = df_data.assign(_yabs=df_data[comparison].abs())
        return ordered.sort_values(by='_yabs', ascending=False).head(annot_abs)
    if annot_top is not None:
        sign = np.sign(annot_top)
        ordered = df_data.sort_values(by=comparison, ascending=(sign == -1))
        return ordered.head(abs(annot_top))
    return df_data.iloc[0:0]


def annotate_points(ax, rows, x_column, y_column, annot_label, offset=0.0,
                    **text_kws):
    """Write a label beside each of the given rows."""
    for _, row in rows.iterrows():
        label = row.get(annot_label)
        if label is None or (isinstance(label, float) and np.isnan(label)):
            continue
        if str(label) == 'nan':
            continue
        ax.text(row[x_column], row[y_column] - offset, label, **text_kws)


def style_axes(ax, xlab=None, ylab=None, title=None,
               xlim_kws=None, ylim_kws=None, hide_spines=('top', 'right')):
    """Labels, limits and spine visibility, applied the same way everywhere."""
    if title is not None:
        ax.set_title(title)
    if xlab is not None:
        ax.set_xlabel(xlab)
    if ylab is not None:
        ax.set_ylabel(ylab)
    if xlim_kws:
        ax.set_xlim(**xlim_kws)
    if ylim_kws:
        ax.set_ylim(**ylim_kws)
    for spine in hide_spines:
        ax.spines[spine].set_visible(False)


def save_and_show(fig, savefig=True, show=True, out_dir='', out_name='plot',
                  out_type='png', dpi=300, tight=True):
    """
    Write the figure out and/or display it.

    Creates out_dir if it does not exist, which previously raised for anyone
    passing a directory they had not made by hand.
    """
    import matplotlib.pyplot as plt

    if tight:
        try:
            fig.tight_layout()
        except Exception:
            pass

    if savefig:
        outpath = Path(out_dir)
        if str(outpath) not in ('', '.'):
            outpath.mkdir(parents=True, exist_ok=True)
        fig.savefig(outpath / f'{out_name}.{out_type}', format=out_type, dpi=dpi)

    if show:
        plt.show()
    else:
        plt.close(fig)


def save_and_show_interactive(fig, savefig=True, show=True, out_dir='',
                              out_name='plot'):
    """The plotly equivalent of save_and_show."""
    if savefig:
        outpath = Path(out_dir)
        if str(outpath) not in ('', '.'):
            outpath.mkdir(parents=True, exist_ok=True)
        fig.write_html(str(outpath / f'{out_name}.html'))
    if show:
        fig.show()
