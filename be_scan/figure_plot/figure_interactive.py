"""
Author: Calvin XiaoYang Hu

{Description: plotly counterparts to every figure_plot figure function}

Each function here mirrors its static twin: same positional arguments, same
AxisLabelOpts / style / OutputOpts objects, same defaults. Switching from a
static figure to its interactive version is a one-word change at the call
site, and the two are styled from the same options, so they agree.

    scatterplot_figure(df, ['score'], 'pos', output=out)
    interactive_scatterplot_figure(df, ['score'], 'pos', output=out)

The HTML is written beside the static file, sharing its stem:

    figures/screen.svg      <- static
    figures/screen.html     <- interactive

Both honour output.save and output.show. Every function also takes `scale`,
which multiplies pixel size and font size together so print-sized figures
stay legible in a browser without changing their proportions.
"""

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from be_scan.figure_plot.figure_classes import *
from be_scan.figure_plot.figure_interactive_style import (
    add_domains, add_hline, apply_axis, base_layout, figure_size,
    finish_interactive, hover_text, line_from_kws, marker_from_kws,
    marker_size, opacity_from_kws, resolve_colors, to_css,
)
from be_scan.plot._annotating_ import calc_neg_ctrls, norm_to_intergenic_ctrls


# --------------------------------------------------------------------------- #
# scatterplot
# --------------------------------------------------------------------------- #


def interactive_scatterplot_figure(
    df_data,
    y_col_list,
    x_col,
    *,
    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    domain: DomainOpts = DOMAIN,
    legend: LegendOpts = LEGEND,
    output: OutputOpts = OUTPUT,
    neg_ctrl: NegCtrlOpts = NEGCONTROL,
    scale: float = 2.0,
    hover_cols=None,
):
    """Interactive twin of scatterplot_figure."""
    df_data = load_frame(df_data)

    needed = set(y_col_list + [x_col])
    missing = needed.difference(df_data.columns)
    if missing:
        raise ValueError(f"Missing column(s): {sorted(missing)}")

    # NEGATIVE CONTROL NORMALISATION, AS IN THE STATIC VERSION #
    if neg_ctrl.adjust:
        assert neg_ctrl.neg_ctrl_col in df_data.columns.tolist()
        _, _, avg_dict = calc_neg_ctrls(
            df_data, y_col_list, neg_ctrl.neg_ctrl_col, neg_ctrl.neg_ctrl_conditions)
        df_data = norm_to_intergenic_ctrls(df_data, y_col_list, avg_dict)

    if axis.xwindow:
        df_data = df_data[df_data[x_col].between(*axis.xwindow)]

    n = len(y_col_list)
    fig = make_subplots(rows=n, cols=1, subplot_titles=y_col_list if n > 1 else None,
                        shared_xaxes=True, vertical_spacing=0.08)

    colors = {}
    if style.hue_col:
        colors = resolve_colors(df_data[style.hue_col].dropna().unique(), style)

    for i, ycol in enumerate(y_col_list, start=1):
        add_domains(fig, domain, df_data, row=i, col=1)

        # the static figure splits points into a faded and an opaque layer at
        # an absolute threshold; the same split is reproduced here
        opaque_mask = df_data[ycol].abs() >= style.transparency_threshold
        layers = [(df_data[~opaque_mask], style.transparent_kws, 'below threshold'),
                  (df_data[opaque_mask], style.opaque_kws, 'above threshold')]

        for layer_df, layer_kws, layer_name in layers:
            if layer_df.empty:
                continue
            groups = ([(v, layer_df[layer_df[style.hue_col] == v])
                       for v in sorted(layer_df[style.hue_col].dropna().unique())]
                      if style.hue_col else [(layer_name, layer_df)])
            for value, sub in groups:
                if sub.empty:
                    continue
                marker = marker_from_kws(
                    layer_kws, scale,
                    color=to_css(colors.get(value)) if style.hue_col else 'black')
                text = hover_text(sub, hover_cols or [style.hue_col, style.marker_col])
                fig.add_trace(go.Scattergl(
                    x=sub[x_col], y=sub[ycol], mode='markers',
                    name=str(value), legendgroup=str(value),
                    showlegend=(i == 1 and layer_name == 'above threshold'
                                or (i == 1 and not style.hue_col and layer_name == 'above threshold')),
                    marker=marker, opacity=opacity_from_kws(layer_kws),
                    hovertext=text,
                    hovertemplate=f'{x_col}: %{{x}}<br>{ycol}: %{{y:.3f}}'
                                  + ('<br>%{hovertext}' if text is not None else '')
                                  + '<extra></extra>',
                ), row=i, col=1)

        if neg_ctrl.lines:
            add_hline(fig, 2 * neg_ctrl.sd, neg_ctrl.ctrl_line_kws, scale, row=i, col=1)
            add_hline(fig, -2 * neg_ctrl.sd, neg_ctrl.ctrl_line_kws, scale, row=i, col=1)

        apply_axis(fig, axis, scale, row=i, col=1, transparent=output.transparent)
        if axis.ylabel is None:
            fig.update_yaxes(title_text=ycol, row=i, col=1)

    fig.update_layout(**base_layout(output, scale, n_rows=n),
                      showlegend=bool(style.hue_col) and legend.loc is not None
                      or bool(style.hue_col))
    return finish_interactive(fig, output)


# --------------------------------------------------------------------------- #
# correlation scatterplot
# --------------------------------------------------------------------------- #


def interactive_correlation_scatterplot_figure(
    df_data,
    x_column,
    y_column,
    legend_order,

    hue_col='',
    palette={},
    highlight_col="mut_type",
    highlight_categories=["Splice", "Intron", "Nonsense"],

    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    hover_cols=None,
    ):
    """Interactive twin of correlation_scatterplot_figure."""
    df_data = load_frame(df_data)

    thr = style.transparency_threshold
    faded = (df_data[x_column].between(-thr, thr)
             & df_data[y_column].between(-thr, thr))

    fig = go.Figure()

    # CUTOFF LINES, MATCHING THE STATIC CROSSHAIRS #
    cutoff_line = line_from_kws(style.cutoff_line_kws, scale)
    for cutoff in (-thr, thr):
        fig.add_hline(y=cutoff, line=cutoff_line)
        fig.add_vline(x=cutoff, line=cutoff_line)

    highlighted = (df_data[highlight_col].isin(highlight_categories)
                   if highlight_col in df_data.columns
                   else pd.Series(False, index=df_data.index))

    layers = [
        (~highlighted & ~faded, style.opaque_kws, True),
        (~highlighted & faded, style.transparent_kws, False),
        (highlighted, style.highlight_kws, True),
    ]

    values = (sorted(df_data[hue_col].dropna().unique()) if hue_col else [None])
    colors = resolve_colors(values, style, palette) if hue_col else {}

    for mask, kws, show in layers:
        sub_all = df_data[mask]
        if sub_all.empty:
            continue
        groups = ([(v, sub_all[sub_all[hue_col] == v]) for v in values]
                  if hue_col else [(None, sub_all)])
        for value, sub in groups:
            if sub.empty:
                continue
            text = hover_text(sub, hover_cols or [hue_col, highlight_col])
            fig.add_trace(go.Scattergl(
                x=sub[x_column], y=sub[y_column], mode='markers',
                name=str(value) if value is not None else y_column,
                legendgroup=str(value), showlegend=show and hue_col != '',
                marker=marker_from_kws(kws, scale,
                                       color=to_css(colors.get(value))),
                opacity=opacity_from_kws(kws),
                hovertext=text,
                hovertemplate=f'{x_column}: %{{x:.3f}}<br>{y_column}: %{{y:.3f}}'
                              + ('<br>%{hovertext}' if text is not None else '')
                              + '<extra></extra>',
            ))

    # LINE OF BEST FIT AND PEARSON, AS IN THE STATIC VERSION #
    x_vals = df_data[x_column].to_numpy(dtype=float)
    y_vals = df_data[y_column].to_numpy(dtype=float)
    ok = ~np.isnan(x_vals) & ~np.isnan(y_vals)
    if ok.sum() > 1:
        slope, intercept = np.polyfit(x_vals[ok], y_vals[ok], 1)
        x_fit = np.linspace(np.min(x_vals[ok]), np.max(x_vals[ok]), 100)
        fig.add_trace(go.Scatter(
            x=x_fit, y=slope * x_fit + intercept, mode='lines',
            line=line_from_kws(style.line_kws, scale),
            name='best fit', showlegend=False, hoverinfo='skip'))
        r_val, p_val = pearsonr(x_vals[ok], y_vals[ok])
        fig.add_annotation(
            xref='paper', yref='paper', x=0.02, y=0.98, showarrow=False,
            align='left', text=f'r = {r_val:.3f}<br>p = {p_val:.2e}')

    fig.update_layout(**base_layout(output, scale))
    apply_axis(fig, axis, scale, transparent=output.transparent)
    if axis.xlabel is None:
        fig.update_xaxes(title_text=x_column)
    if axis.ylabel is None:
        fig.update_yaxes(title_text=y_column)
    return finish_interactive(fig, output)


# --------------------------------------------------------------------------- #
# boxplot
# --------------------------------------------------------------------------- #


def interactive_boxplot_figure(
    df,
    x_col, y_col,
    palette,
    lab_order=None,

    axis: AxisLabelOpts = AXIS,
    style: BoxStyleOpts = BOXPLOTSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """Interactive twin of boxplot_figure."""
    df_data = load_frame(df)
    order = lab_order if lab_order is not None else sorted(df_data[x_col].dropna().unique())

    fig = go.Figure()
    add_hline(fig, 0, style.axhline_kws, scale)

    box_kws = style.box_kws or {}
    for value in order:
        sub = df_data[df_data[x_col] == value]
        if sub.empty:
            continue
        fig.add_trace(go.Box(
            y=sub[y_col], name=str(value),
            marker_color=to_css((palette or {}).get(value)),
            width=box_kws.get('width', 0.4),
            boxpoints='outliers' if box_kws.get('showfliers', False) else False,
            line=dict(width=style.boxprops.get('linewidth', 0.5) * scale,
                      color=to_css(style.boxprops.get('edgecolor', 'black'))),
            fillcolor=to_css((palette or {}).get(value)),
            hovertemplate=f'{x_col}: {value}<br>{y_col}: %{{y:.3f}}<extra></extra>',
        ))

    fig.update_layout(**base_layout(output, scale), showlegend=False,
                      boxmode='group')
    apply_axis(fig, axis, scale, transparent=output.transparent)
    fig.update_xaxes(tickangle=-style.xlabel_rotation)
    if axis.ylabel is None:
        fig.update_yaxes(title_text=y_col)
    return finish_interactive(fig, output)


# --------------------------------------------------------------------------- #
# jitter box + kde
# --------------------------------------------------------------------------- #


def interactive_jitterbox_kdeplot_figure(
    data,
    xcol, ycol,
    huecol, color_dict,
    jitterbox=True, kde=True, threshold=None,

    axis: AxisLabelOpts = AXIS,
    style: JitterStyleOpts = JITTERPLOTSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """
    Interactive twin of jitterbox_kdeplot_figure.

    Returns a dict keyed like the static version, so the two can be swapped
    without changing how the result is unpacked.
    """
    data = load_frame(data)
    figures = {}

    if jitterbox:
        fig = go.Figure()
        for value in sorted(data[xcol].dropna().unique()):
            sub = data[data[xcol] == value]
            fig.add_trace(go.Box(
                y=sub[ycol], name=str(value),
                boxpoints=False, fillcolor='rgba(0,0,0,0)',
                line=dict(color=to_css(style.whiskerprops.get('color', 'gray')),
                          width=style.boxprops.get('linewidth', 0.5) * scale),
                hoverinfo='skip', showlegend=False,
            ))
        # the static version splits points at `threshold` into faded and solid
        for hue_value in sorted(data[huecol].dropna().unique()):
            sub = data[data[huecol] == hue_value]
            below = sub[sub[ycol] <= threshold] if threshold is not None else sub
            above = sub[sub[ycol] > threshold] if threshold is not None else sub.iloc[0:0]
            for part, kws, show in ((below, style.stripplot_kws, True),
                                    (above, style.stripplot_hits_kws, False)):
                if part.empty:
                    continue
                fig.add_trace(go.Box(
                    y=part[ycol], x=part[xcol].astype(str),
                    name=str(hue_value), legendgroup=str(hue_value),
                    showlegend=show, boxpoints='all', pointpos=0, jitter=0.5,
                    fillcolor='rgba(0,0,0,0)', line=dict(width=0),
                    marker=dict(color=to_css((color_dict or {}).get(hue_value)),
                                size=marker_size(kws.get('size', 2) ** 2, scale),
                                opacity=kws.get('alpha', 0.75)),
                    hovertemplate=f'{ycol}: %{{y:.3f}}<extra>{hue_value}</extra>',
                ))
        fig.update_layout(**base_layout(output, scale), boxmode='overlay')
        apply_axis(fig, axis, scale, transparent=output.transparent)
        figures['jitterbox'] = finish_interactive(fig, output, suffix='_Boxplot')

    if kde:
        from scipy.stats import gaussian_kde
        fig = go.Figure()
        quartiles = sorted(data[xcol].dropna().unique())[::-1]
        colors = list(style.colors)[::-1]
        for q, colour in zip(quartiles, colors):
            values = data.loc[data[xcol] == q, ycol].dropna().to_numpy(dtype=float)
            if values.size < 2:
                continue
            grid = np.linspace(values.min(), values.max(), 200)
            density = gaussian_kde(values)(grid)
            fig.add_trace(go.Scatter(
                x=density, y=grid, mode='lines', name=f'Q{q}',
                line=dict(color=to_css(colour),
                          width=style.kdeplot_kws.get('linewidth', 0.5) * scale),
                hovertemplate=f'{ycol}: %{{y:.3f}}<br>density: %{{x:.3f}}'
                              f'<extra>Q{q}</extra>',
            ))
        fig.update_layout(**base_layout(output, scale))
        apply_axis(fig, axis, scale, transparent=output.transparent)
        fig.update_xaxes(title_text='Density')
        figures['kde'] = finish_interactive(fig, output, suffix='_Density')

    return figures


# --------------------------------------------------------------------------- #
# lollipop
# --------------------------------------------------------------------------- #


def _lollipop_traces(fig, df, xcol, ycol, threshold, category_col, colors,
                     style, scale, row=None, col=None, hover_cols=None):
    """Stems and markers for one half of a lollipop plot."""
    if df is None or len(df) == 0:
        return
    df = df.copy()
    hits = df[df[ycol].abs() >= threshold] if threshold is not None else df

    xs, ys = [], []
    for _, r in hits.iterrows():
        xs += [r[xcol], r[xcol], None]
        ys += [0, r[ycol], None]
    stem = go.Scatter(x=xs, y=ys, mode='lines', showlegend=False,
                      hoverinfo='skip',
                      line=dict(color=to_css(style.stem_kws.get('color', 'black')),
                                width=style.stem_kws.get('linewidth', 0.4) * scale))
    fig.add_trace(stem, row=row, col=col) if row else fig.add_trace(stem)

    groups = (hits.groupby(category_col) if category_col in hits.columns
              else [(None, hits)])
    for value, sub in groups:
        text = hover_text(sub, hover_cols or [category_col])
        marker = go.Scatter(
            x=sub[xcol], y=sub[ycol], mode='markers', name=str(value),
            legendgroup=str(value),
            marker=dict(size=style.marker_kws.get('markersize', 3) * scale,
                        color=to_css((colors or {}).get(value)),
                        line=dict(color=to_css(style.marker_kws.get('markeredgecolor', 'black')),
                                  width=style.marker_kws.get('markeredgewidth', 0.4) * scale)),
            hovertext=text,
            hovertemplate=f'{xcol}: %{{x}}<br>{ycol}: %{{y:.3f}}'
                          + ('<br>%{hovertext}' if text is not None else '')
                          + '<extra></extra>',
        )
        fig.add_trace(marker, row=row, col=col) if row else fig.add_trace(marker)


def interactive_lollipop_figure(
    df_pos, pos_xcol, pos_ycol, pos_threshold, pos_category_col, pos_colors,
    df_neg, neg_xcol, neg_ycol, neg_threshold, neg_category_col, neg_colors,

    axis: AxisLabelOpts = AXIS,
    domain: DomainOpts = DOMAIN,
    style: LollipopStyleOpts = LOLLIPOPSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    hover_cols=None,
    ):
    """Interactive twin of lollipop_figure."""
    fig = go.Figure()
    add_domains(fig, domain, load_frame(df_pos) if df_pos is not None else None)
    add_hline(fig, 0, style.base_kws, scale)

    _lollipop_traces(fig, load_frame(df_pos), pos_xcol, pos_ycol, pos_threshold,
                     pos_category_col, pos_colors, style, scale, hover_cols=hover_cols)
    _lollipop_traces(fig, load_frame(df_neg), neg_xcol, neg_ycol, neg_threshold,
                     neg_category_col, neg_colors, style, scale, hover_cols=hover_cols)

    for thr, sign in ((pos_threshold, 1), (neg_threshold, -1)):
        if thr is not None:
            add_hline(fig, sign * abs(thr), style.line_kws, scale)

    fig.update_layout(**base_layout(output, scale))
    apply_axis(fig, axis, scale, transparent=output.transparent)
    return finish_interactive(fig, output)


def interactive_lollipop_split_figure(
    df_pos, pos_xcol, pos_ycol, pos_threshold, pos_category_col, pos_colors,
    ylim_top, yticks_top,
    df_neg, neg_xcol, neg_ycol, neg_threshold, neg_category_col, neg_colors,
    ylim_bot, yticks_bot,
    height_ratios=[1, 0.001, 1],

    axis: AxisLabelOpts = AXIS,
    domain: DomainOpts = DOMAIN,
    style: LollipopStyleOpts = LOLLIPOPSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    hover_cols=None,
    ):
    """Interactive twin of lollipop_split_figure, with a broken y axis."""
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                        vertical_spacing=0.03,
                        row_heights=[height_ratios[0], height_ratios[-1]])

    add_domains(fig, domain, load_frame(df_pos) if df_pos is not None else None,
                row=1, col=1)
    add_domains(fig, domain, load_frame(df_neg) if df_neg is not None else None,
                row=2, col=1)

    _lollipop_traces(fig, load_frame(df_pos), pos_xcol, pos_ycol, pos_threshold,
                     pos_category_col, pos_colors, style, scale,
                     row=1, col=1, hover_cols=hover_cols)
    _lollipop_traces(fig, load_frame(df_neg), neg_xcol, neg_ycol, neg_threshold,
                     neg_category_col, neg_colors, style, scale,
                     row=2, col=1, hover_cols=hover_cols)

    fig.update_layout(**base_layout(output, scale, n_rows=2))
    apply_axis(fig, axis, scale, row=1, col=1, transparent=output.transparent)
    apply_axis(fig, axis, scale, row=2, col=1, transparent=output.transparent)
    if ylim_top is not None:
        fig.update_yaxes(range=list(ylim_top), row=1, col=1)
    if ylim_bot is not None:
        fig.update_yaxes(range=list(ylim_bot), row=2, col=1)
    if yticks_top is not None:
        fig.update_yaxes(tickmode='array', tickvals=list(yticks_top), row=1, col=1)
    if yticks_bot is not None:
        fig.update_yaxes(tickmode='array', tickvals=list(yticks_bot), row=2, col=1)
    return finish_interactive(fig, output)


# --------------------------------------------------------------------------- #
# clustered heatmap
# --------------------------------------------------------------------------- #


def interactive_kmeans_clustered_heatmap(
    df, k,
    clusters,
    genes,
    screen_order=None,
    gene_palette={},
    gene_list=None,
    cluster_order=None,
    highlight_cols=None,
    order_col=None,

    axis: AxisLabelOpts = AXIS,
    style: ClusteredHeatmapStyleOpts = CLUSTEREDHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """Interactive twin of kmeans_clustered_heatmap."""
    df_data = load_frame(df)
    order = screen_order if screen_order is not None else list(df_data.columns)
    matrix = df_data[order]

    labels = pd.Series(list(clusters), index=df_data.index, name='cluster')
    sort_keys = [labels]
    if order_col is not None and order_col in df_data.columns:
        sort_keys.append(df_data[order_col])
    ordering = pd.concat(sort_keys, axis=1).sort_values(
        list(pd.concat(sort_keys, axis=1).columns)).index
    matrix = matrix.loc[ordering]

    heat = style.heatmap_kws or {}
    fig = go.Figure(go.Heatmap(
        z=matrix.to_numpy(dtype=float),
        x=list(matrix.columns),
        y=[str(i) for i in matrix.index],
        colorscale='RdBu_r' if heat.get('cmap') == 'bwr' else 'RdBu_r',
        zmid=heat.get('center', 0),
        zmin=heat.get('vmin'), zmax=heat.get('vmax'),
        showscale=heat.get('cbar', False) is not False,
        hovertemplate='%{x}<br>guide %{y}<br>score %{z:.3f}<extra></extra>',
    ))

    # cluster boundaries, matching the static border lines
    boundaries = np.flatnonzero(np.diff(labels.loc[ordering].to_numpy())) + 0.5
    for b in boundaries:
        fig.add_hline(y=b, line=dict(color=to_css(style.line_kws.get('color', 'black')),
                                     width=style.line_kws.get('linewidth', 0.75) * scale))

    fig.update_layout(**base_layout(output, scale))
    apply_axis(fig, axis, scale, transparent=output.transparent)
    fig.update_yaxes(showticklabels=False, ticks='')
    return finish_interactive(fig, output)


# --------------------------------------------------------------------------- #
# loess
# --------------------------------------------------------------------------- #


def interactive_loess_figure(
    results: dict,
    gene: str = '',
    max_cap: int = 4,
    pval_thresh: float = 0.05,
    *,
    axis_pval: AxisLabelOpts = AXIS,
    axis_regr: AxisLabelOpts = AXIS,
    style: LoessStyleOpts = LOESSSTYLE,
    domain: DomainOpts = DOMAIN,
    output_pval: OutputOpts = OUTPUT,
    output_regr: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """
    Interactive twin of loess_figure.

    Returns (pval_figure, regression_figure), matching the pair the static
    version produces.
    """
    threshold = -np.log10(pval_thresh)

    # ---- significance track ---- #
    fig_pval = go.Figure()
    add_domains(fig_pval, domain)
    for name, res in results.items():
        x = np.asarray(res['x_obs' if 'x_obs' in res else 'x'])
        y = np.asarray(res['y_pval' if 'y_pval' in res else 'pval'])
        y = np.minimum(y, max_cap)
        above = y >= threshold
        for mask, colour, label in ((~above, style.below_color, 'n.s.'),
                                    (above, style.above_color, 'significant')):
            if not mask.any():
                continue
            fig_pval.add_trace(go.Scattergl(
                x=x[mask], y=y[mask], mode='markers', name=f'{name} {label}',
                marker=dict(size=style.pval_plot_kws.get('markersize', 1) * scale * 2,
                            color=to_css(colour)),
                opacity=style.pval_plot_kws.get('alpha', 1),
                hovertemplate='position %{x}<br>-log10(p) %{y:.2f}<extra></extra>',
            ))
    fig_pval.add_hline(y=threshold, line=line_from_kws(style.line_kws, scale))
    fig_pval.update_layout(**base_layout(output_pval, scale))
    apply_axis(fig_pval, axis_pval, scale, transparent=output_pval.transparent)
    if axis_pval.ylabel is None:
        fig_pval.update_yaxes(title_text='-log10(p)')
    finish_interactive(fig_pval, output_pval, suffix='_pval')

    # ---- regression track ---- #
    fig_regr = go.Figure()
    add_domains(fig_regr, domain)
    for name, res in results.items():
        x = np.asarray(res['x_obs' if 'x_obs' in res else 'x'])
        y = np.asarray(res['y_loess' if 'y_loess' in res else 'loess'])
        fig_regr.add_trace(go.Scatter(
            x=x, y=y, mode='lines', name=str(name),
            line=dict(color=to_css(style.regr_plot_kws.get('color', 'gray')),
                      width=scale),
            hovertemplate='position %{x}<br>smoothed %{y:.3f}<extra></extra>',
        ))
    fig_regr.update_layout(**base_layout(output_regr, scale))
    apply_axis(fig_regr, axis_regr, scale, transparent=output_regr.transparent)
    finish_interactive(fig_regr, output_regr, suffix='_regr')

    return fig_pval, fig_regr


# --------------------------------------------------------------------------- #
# PWES
# --------------------------------------------------------------------------- #


def interactive_pwes_histogram(
    df_clus,
    color_list,

    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """Interactive twin of pwes_histogram."""
    df_data = load_frame(df_clus)
    counts = df_data['cluster'].value_counts().sort_index()
    colors = [color_list[(i - 1) % len(color_list)] for i in counts.index]

    fig = go.Figure(go.Bar(
        x=counts.index.astype(str), y=counts.to_numpy(),
        marker=dict(color=[to_css(c) for c in colors], line=dict(color='black', width=0.5 * scale)),
        hovertemplate='cluster %{x}<br>%{y} guides<extra></extra>'))
    fig.update_layout(**base_layout(output, scale), showlegend=False)
    apply_axis(fig, axis, scale, transparent=output.transparent)
    fig.update_xaxes(title_text=axis.xlabel or 'Cluster Number')
    fig.update_yaxes(title_text=axis.ylabel or 'Count')
    return finish_interactive(fig, output, suffix='-cluster_histogram')


def interactive_pwes_boxplots(
    df_clus,
    x_col,
    color_list,

    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """Interactive twin of pwes_boxplots."""
    df_data = load_frame(df_clus)
    fig = go.Figure()
    for i, (cluster, sub) in enumerate(df_data.groupby('cluster')):
        fig.add_trace(go.Box(
            y=sub[x_col], name=str(cluster), boxpoints='all', jitter=0.4,
            pointpos=0, fillcolor='rgba(0,0,0,0)',
            line=dict(color='black', width=0.5 * scale),
            marker=dict(color=to_css(color_list[i % len(color_list)]),
                        size=2 * scale, opacity=0.5),
            hovertemplate=f'cluster {cluster}<br>{x_col}: %{{y:.3f}}<extra></extra>'))
    fig.update_layout(**base_layout(output, scale), showlegend=False)
    apply_axis(fig, axis, scale, transparent=output.transparent)
    fig.update_xaxes(title_text=axis.xlabel or 'Cluster Number')
    fig.update_yaxes(title_text=axis.ylabel or f'{x_col} Values')
    return finish_interactive(fig, output, suffix='-cluster_boxplots')


def interactive_pwes_scatterplots(
    df_clus,
    x_col, scores_col,
    color_list,
    max_cluster=None,

    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    domain: DomainOpts = DOMAIN,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """Interactive twin of pwes_scatterplots."""
    df_data = load_frame(df_clus)
    fig = go.Figure()
    add_domains(fig, domain, df_data)
    for i, (cluster, sub) in enumerate(df_data.groupby('cluster')):
        if max_cluster is not None and cluster > max_cluster:
            continue
        fig.add_trace(go.Scattergl(
            x=sub[x_col], y=sub[scores_col], mode='markers',
            name=f'cluster {cluster}',
            marker=marker_from_kws(style.opaque_kws, scale,
                                   color=to_css(color_list[i % len(color_list)])),
            opacity=opacity_from_kws(style.opaque_kws),
            hovertemplate=f'{x_col}: %{{x}}<br>{scores_col}: %{{y:.3f}}'
                          f'<extra>cluster {cluster}</extra>'))
    fig.update_layout(**base_layout(output, scale))
    apply_axis(fig, axis, scale, transparent=output.transparent)
    return finish_interactive(fig, output, suffix='-cluster_scatterplots')


def interactive_pwes_heatmap(
    df_scaled,
    mask_on=True,
    blackout=['X'],

    axis: AxisLabelOpts = AXIS,
    style: PWESHeatmapStyleOpts = PWESHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """Interactive twin of pwes_heatmap."""
    matrix = load_frame(df_scaled)
    heat = style.heatmap_kws or {}
    fig = go.Figure(go.Heatmap(
        z=matrix.to_numpy(dtype=float),
        x=[str(c) for c in matrix.columns],
        y=[str(i) for i in matrix.index],
        colorscale='RdBu_r',
        zmin=heat.get('vmin', -1), zmax=heat.get('vmax', 1), zmid=0,
        showscale=bool(heat.get('cbar', False)),
        hovertemplate='%{y} vs %{x}<br>PWES %{z:.3f}<extra></extra>'))
    fig.update_layout(**base_layout(output, scale))
    apply_axis(fig, axis, scale, transparent=output.transparent)
    fig.update_xaxes(showticklabels=bool(heat.get('xticklabels', False)))
    fig.update_yaxes(showticklabels=bool(heat.get('yticklabels', False)),
                     autorange='reversed')
    return finish_interactive(fig, output, suffix='-PWES_heatmap')


def interactive_pwes_clustermap(
    df_scaled, link, df_clusters,
    color_list=[],

    axis: AxisLabelOpts = AXIS,
    style: PWESHeatmapStyleOpts = PWESHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """
    Interactive twin of pwes_clustermap.

    Plotly has no clustermap, so the rows and columns are reordered by the
    same linkage the static version draws its dendrogram from, and cluster
    boundaries are marked instead.
    """
    from scipy.cluster.hierarchy import leaves_list

    matrix = load_frame(df_scaled)
    try:
        order = leaves_list(link)
    except Exception:
        order = np.arange(len(matrix))
    ordered = matrix.iloc[order, order] if matrix.shape[0] == matrix.shape[1] \
        else matrix.iloc[order]

    heat = style.clustered_heatmap_kws or {}
    fig = go.Figure(go.Heatmap(
        z=ordered.to_numpy(dtype=float),
        colorscale='RdBu_r',
        zmin=heat.get('vmin', -1), zmax=heat.get('vmax', 1),
        zmid=heat.get('center', 0), showscale=True,
        hovertemplate='PWES %{z:.3f}<extra></extra>'))
    fig.update_layout(**base_layout(output, scale))
    apply_axis(fig, axis, scale, transparent=output.transparent)
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False, autorange='reversed')
    return finish_interactive(fig, output, suffix='-cluster_clustermap')


def interactive_pwes_heatmap_colorbar(
    orientation='vertical',

    axis: AxisLabelOpts = AXIS,
    style: PWESHeatmapStyleOpts = PWESHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    scale: float = 2.0,
    ):
    """Interactive twin of pwes_heatmap_colorbar: the colour scale alone."""
    norm = style.colorbar_normalize_kws or {}
    vmin, vmax = norm.get('vmin', -1), norm.get('vmax', 1)
    grid = np.linspace(vmin, vmax, 256)
    z = grid.reshape(-1, 1) if orientation.startswith('v') else grid.reshape(1, -1)

    fig = go.Figure(go.Heatmap(
        z=z, colorscale='RdBu_r', zmin=vmin, zmax=vmax, zmid=0,
        showscale=False, hovertemplate='%{z:.2f}<extra></extra>'))
    fig.update_layout(**base_layout(output, scale))
    apply_axis(fig, axis, scale, transparent=output.transparent)
    fig.update_xaxes(showticklabels=orientation.startswith('h'))
    fig.update_yaxes(showticklabels=orientation.startswith('v'))
    return finish_interactive(fig, output, suffix='-colorbar')
