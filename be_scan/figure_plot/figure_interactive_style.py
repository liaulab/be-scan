"""
Author: Calvin XiaoYang Hu

{Description: shared plotly styling so interactive figures match the static ones}

Every figure_plot function has an interactive counterpart. The pair is driven
by the same AxisLabelOpts / style / OutputOpts objects, and the plotly layout
below is derived from the same matplotlib rcParams the static figures use, so
the two renderings agree rather than merely resembling each other.

What is carried across:

  font          Arial, falling back to sans-serif, at the rcParams size
  spines        left and bottom only, at rcParams axes.linewidth
  ticks         outside, at rcParams tick width and length
  background    transparent when OutputOpts.transparent is set
  gridlines     off, as in the static figures
  size          OutputOpts.figsize, which is in inches, converted to pixels
  colours       resolved through the same palette logic as the static plots

One thing deliberately does not carry across: absolute size. The static
figures are sized for print, typically 4 x 2.2 inches with 6 pt text, which
is illegible in a browser. Every interactive function therefore takes a
`scale` argument that multiplies pixel dimensions and font sizes by the same
factor, so proportions are preserved exactly while the result is readable.
scale=1.0 reproduces the print geometry.
"""

import math

import matplotlib.colors as mcolors
import numpy as np
import plotly.graph_objects as go

from be_scan.figure_plot.figure_classes import (
    AXIS, OUTPUT, mpl, resolve_colors, resolve_domains,
)

# matplotlib reports figsize in inches; this is the conversion plotly needs
PX_PER_INCH = 96.0

# how many points of font the static figures use, read from the live rcParams
# so restyling the static plots restyles the interactive ones too
BASE_FONT = float(mpl.rcParams.get('font.size', 6))
BASE_AXES_LW = float(mpl.rcParams.get('axes.linewidth', 0.5))
BASE_TICK_W = float(mpl.rcParams.get('xtick.major.width', 0.5))
BASE_TICK_LEN = float(mpl.rcParams.get('xtick.major.size', 2))
FONT_FAMILY = 'Arial, Helvetica, sans-serif'

# How plotly.js is delivered in the written HTML. True embeds it, giving a
# self-contained ~3.7 MB file that opens offline; 'cdn' loads it from the
# network instead and shrinks each file to tens of kilobytes.
INCLUDE_PLOTLYJS = True


def figure_size(output, scale=2.0):
    """Pixel width and height matching OutputOpts.figsize, in inches."""
    width, height = output.figsize
    return width * PX_PER_INCH * scale, height * PX_PER_INCH * scale


def axis_style(scale=1.0, transparent=True):
    """The per-axis settings shared by every interactive figure."""
    return dict(
        showline=True,
        linewidth=BASE_AXES_LW * scale,
        linecolor='black',
        mirror=False,               # no top/right spine, as in the static plots
        showgrid=False,
        zeroline=False,
        ticks='outside',
        tickwidth=BASE_TICK_W * scale,
        ticklen=BASE_TICK_LEN * scale,
        tickfont=dict(family=FONT_FAMILY, size=BASE_FONT * scale, color='black'),
        title=dict(font=dict(family=FONT_FAMILY, size=BASE_FONT * scale,
                             color='black')),
    )


def base_layout(output, scale=2.0, n_rows=1):
    """The layout every interactive figure starts from."""
    width, height = figure_size(output, scale)
    background = 'rgba(0,0,0,0)' if output.transparent else 'white'
    return dict(
        width=width,
        height=height * n_rows,
        paper_bgcolor=background,
        plot_bgcolor=background,
        font=dict(family=FONT_FAMILY, size=BASE_FONT * scale, color='black'),
        margin=dict(l=50 * scale, r=20 * scale, t=30 * scale, b=40 * scale),
        title=dict(font=dict(family=FONT_FAMILY, size=BASE_FONT * scale)),
        legend=dict(font=dict(family=FONT_FAMILY, size=BASE_FONT * scale),
                    bgcolor='rgba(0,0,0,0)'),
        hoverlabel=dict(font=dict(family=FONT_FAMILY, size=BASE_FONT * scale)),
    )


def apply_axis(fig, axis, scale=2.0, row=None, col=None, transparent=True):
    """
    Push an AxisLabelOpts onto a figure, the same way the static plots do.

    Only the fields that are set are applied, matching apply_axis_limits on
    the static side.
    """
    style = axis_style(scale, transparent)

    xaxis = dict(style)
    yaxis = dict(style)
    if axis.xlabel is not None:
        xaxis['title']['text'] = axis.xlabel
    if axis.ylabel is not None:
        yaxis['title']['text'] = axis.ylabel
    if axis.xlim is not None:
        xaxis['range'] = list(axis.xlim)
    if axis.ylim is not None:
        yaxis['range'] = list(axis.ylim)
    if axis.xticks is not None:
        xaxis.update(tickmode='array', tickvals=list(axis.xticks))
    elif axis.protein_len:
        xaxis.update(tickmode='array', tickvals=[1, axis.protein_len])
    if axis.yticks is not None:
        yaxis.update(tickmode='array', tickvals=list(axis.yticks))

    if row is None:
        fig.update_xaxes(**xaxis)
        fig.update_yaxes(**yaxis)
    else:
        fig.update_xaxes(row=row, col=col or 1, **xaxis)
        fig.update_yaxes(row=row, col=col or 1, **yaxis)

    if axis.title:
        fig.update_layout(title=dict(text=axis.title))
    return fig


def marker_size(mpl_size, scale=2.0):
    """
    Convert a matplotlib scatter size to a plotly marker diameter.

    Matplotlib's `s` is an area in points squared, plotly's marker size is a
    diameter, so the conversion is a square root rather than a ratio. Without
    this the interactive points come out wildly larger than the static ones.
    """
    if mpl_size is None:
        return 6 * scale
    return max(math.sqrt(float(mpl_size)) * scale, 1.0)


def marker_from_kws(kws, scale=2.0, color=None):
    """
    Translate a matplotlib scatter kwargs dict into a plotly marker dict.

    The static plots carry their point styling in dataclass fields such as
    ScatterStyleOpts.transparent_kws and .opaque_kws. Reading those same
    dicts here is what keeps the two renderings matched.
    """
    kws = kws or {}
    marker = dict(size=marker_size(kws.get('s'), scale))
    if color is not None:
        marker['color'] = to_css(color)
    edge = kws.get('edgecolor')
    width = kws.get('linewidth', 0)
    if edge and width:
        marker['line'] = dict(color=to_css(edge), width=float(width) * scale)
    else:
        marker['line'] = dict(width=0)
    return marker


def to_css(color, default='black'):
    """
    Render a matplotlib colour as something plotly accepts.

    The style dataclasses use matplotlib conventions, including single-letter
    shorthands: LoessStyleOpts.line_kws and NegCtrlOpts.ctrl_line_kws both
    specify "k" for black. Plotly rejects those, so every colour taken from a
    style object is passed through here.
    """
    if color is None:
        return default
    try:
        return mcolors.to_hex(color)
    except (ValueError, TypeError):
        return default


def line_from_kws(kws, scale=2.0):
    """Translate matplotlib line kwargs into a plotly line dict."""
    kws = kws or {}
    dash_map = {'--': 'dash', ':': 'dot', '-.': 'dashdot', '-': 'solid'}
    style = kws.get('ls', kws.get('linestyle', '-'))
    return dict(
        color=to_css(kws.get('color', kws.get('c'))),
        width=float(kws.get('lw', kws.get('linewidth', 1))) * scale,
        dash=dash_map.get(style, 'solid'),
    )


def opacity_from_kws(kws):
    """The alpha a matplotlib kwargs dict asks for."""
    return float((kws or {}).get('alpha', 1.0))


def add_domains(fig, domain, df=None, row=None, col=None):
    """Shade protein domains, matching the static axvspan bands."""
    for dom in resolve_domains(domain, df):
        kws = dict(
            type='rect', xref='x', yref='paper',
            x0=dom['start'], x1=dom['end'], y0=0, y1=1,
            fillcolor=to_css(dom.get('color', domain.default_color)),
            opacity=dom.get('alpha', domain.default_alpha),
            line_width=0, layer='below',
        )
        if row is None:
            fig.add_shape(**kws)
        else:
            fig.add_shape(row=row, col=col or 1, **kws)
    return fig


def add_hline(fig, y, kws, scale=2.0, row=None, col=None):
    """A horizontal reference line matching the static axhline."""
    line = line_from_kws(kws, scale)
    if row is None:
        fig.add_hline(y=y, line=line)
    else:
        fig.add_hline(y=y, line=line, row=row, col=col or 1)
    return fig


def html_path(output, suffix=''):
    """
    Where the interactive figure is written.

    Deliberately the same stem as the static figure, so a plot and its
    interactive twin sit side by side and differ only in extension.
    """
    from be_scan.figure_plot.figure_classes import figure_path
    return figure_path(output, suffix=suffix, ext='html')


def finish_interactive(fig, output, suffix='', include_plotlyjs=None):
    """
    Write the HTML next to the static figure, show it, and return it.

    Mirrors finish_figure on the static side: honours output.save and
    output.show, and returns the figure either way.

    include_plotlyjs defaults to the module setting INCLUDE_PLOTLYJS. True
    embeds plotly.js in the file, which makes each figure about 3.7 MB but
    self-contained, so it opens with no network and survives being emailed to
    a collaborator. Set be_scan.figure_plot.figure_interactive_style
    .INCLUDE_PLOTLYJS = 'cdn' to load the library from a CDN instead and
    shrink every file to a few tens of kilobytes, at the cost of needing a
    connection to view them.
    """
    if include_plotlyjs is None:
        include_plotlyjs = INCLUDE_PLOTLYJS
    if output.save:
        path = html_path(output, suffix)
        fig.write_html(str(path), include_plotlyjs=include_plotlyjs)
    if output.show:
        fig.show()
    return fig


def hover_text(df, columns):
    """Build a hover string from whichever of `columns` the frame has."""
    present = [c for c in columns if c is not None and c in df.columns]
    if not present:
        return None
    return df[present].astype(str).apply(
        lambda row: '<br>'.join(f'{c}: {v}' for c, v in zip(present, row)), axis=1)
