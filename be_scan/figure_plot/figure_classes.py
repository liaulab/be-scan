# COMMON IMPORTS #
import sys
import os
import gc
import math
import time
import random
import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns
from biopandas.pdb import PandasPdb

# MATPLOTLIB #
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import matplotlib.cm as cm
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle
from matplotlib import gridspec

# DATA CLASS IMPORTS #
from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Optional, Tuple, Union
from pathlib import Path

# IGNORE WARNING #
import logging
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)

# SCIPY AND STATS #
import scipy.cluster as sp_cl
import scipy.interpolate as interp
from scipy.stats import pearsonr
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import fcluster

import statsmodels.stats.multitest as smm
from statsmodels.nonparametric.smoothers_lowess import lowess
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

import pickle

# SETUP MPL SETTINGS #

mpl.rcParams.update({
    # Keep text as text in SVG
    'svg.fonttype': 'none',

    # Font
    'font.family': 'sans-serif',
    'font.sans-serif': 'Arial',

    # Font sizes
    'font.size': 6,
    'axes.labelsize': 6,
    'axes.titlesize': 6,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,

    # Thin axes lines
    'axes.linewidth': 0.5,
    'lines.linewidth': 0.5,

    # Thin tick lines
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 2,
    'ytick.major.size': 2,

    # Smaller padding
    'axes.labelpad' : 2,
    'ytick.major.pad' : 1,
    'xtick.major.pad': 1
})

cm_2_in = 0.393
in_2_cm = 2.54
five_across_size = (4.3/in_2_cm, 2.8/in_2_cm)
three_across_size = (6.35/in_2_cm, 2.8/in_2_cm)

# Dataclass for axes, labels, etc
@dataclass(frozen = True)
class AxisLabelOpts:
    title: Optional[str] = None
    xlabel: Optional[str] = None
    ylabel: Optional[str] = None
    xwindow: Optional[Tuple[float, float]] = None
    xlim: Optional[Tuple[float, float]] = None
    ylim: Optional[Tuple[float, float]] = None
    xticks: Optional[List[int]] = None
    yticks: Optional[List[int]] = None
    protein_len: Optional[int] = None
    linewidth: Optional[float] = 0.5
    tick_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"labelsize":8, "width":0.5})
    space: Optional[float] = 3

@dataclass(frozen = True)
class NegCtrlOpts:
    ###
    adjust: bool = False
    neg_ctrl_col: Optional[str] = 'gene'
    neg_ctrl_conditions: Optional[str] = 'NT Controls'
    lines: bool = False
    sd: Optional[float] = 1
    ctrl_line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "k", "ls": "--", "lw": 1})

# Dataclass for style, colors, etc
@dataclass(frozen = True)
class ScatterStyleOpts:
    palette: Union[str, List[str]] = "colorblind"
    hue_col: Optional[str] = None
    marker_col: Optional[str] = None
    marker_map: Optional[Dict[str, str]] = field(default_factory=dict)
    color_map: Optional[Dict[Any, str]] = None
    transparent_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"edgecolor": "black", "alpha": 0.6, "s": 3, "linewidth": 0.1})
    opaque_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"edgecolor": "black", "alpha": 1, "s": 5, "linewidth": 0.25})
    highlight_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"edgecolor": "black", "alpha": 1, "s": 5, "linewidth": 0.3})
    transparency_threshold: Optional[float] = 2.0
    line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "k", "ls": "--", "lw": 1})
    cutoff_line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "gray", "ls": ":", "lw": 0.5})

@dataclass(frozen = True)
class BoxStyleOpts:
    color_map: Optional[Dict[Any, str]] = None
    xlabel_rotation: Optional[float] = 45
    xlabel_ha: Optional[str] = "right"

    axhline_kws: Dict[str, Any] = field(default_factory=lambda:
        {'linestyle': '--', 'color': 'grey'})
    box_kws: Dict[str, Any] = field(default_factory=lambda:
        {"width":0.4, "showfliers":False, "saturation":1})
    boxprops: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':0.5, 'edgecolor':'black'})
    whiskerprops: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':0.5})
    capprops: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':0.5})
    medianprops: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':1.5, 'color':'black'})

@dataclass(frozen = True)
class JitterStyleOpts:
    linewidth: Optional[float] = 0.5
    colors: Optional[List[str]] = field(default_factory=lambda:
        ['red', 'purple', 'blue', 'cornflowerblue'])
    axhline_kws: Dict[str, Any] = field(default_factory=lambda:
        {'linestyle': '--', 'color': 'grey'})

    box_kws: Dict[str, Any] = field(default_factory=lambda:
        {'width':0.3, 'showcaps':True, 'showfliers':False})
    stripplot_kws: Dict[str, Any] = field(default_factory=lambda:
        {'alpha': 0.75, 'size': 2, 'linewidth': 0, 'zorder':0,
         'jitter':True, 'dodge':False, 'rasterized':True})
    stripplot_hits_kws: Dict[str, Any] = field(default_factory=lambda:
        {'alpha': 1.0, 'size': 3, 'linewidth': 0.25, 'zorder':0,
         'jitter':True, 'dodge':False, 'rasterized':False})
    kdeplot_kws: Dict[str, Any] = field(default_factory=lambda:
        {'linewidth':0.5})

    boxprops: Dict[str, Any] = field(default_factory=lambda:
        {'facecolor': 'none', 'zorder': 1, 'linewidth': 0.5})
    whiskerprops: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'gray', 'zorder': 1, 'linewidth': 0.5})
    capprops: Dict[str, Any] = field(default_factory=lambda:
        {'zorder': 1, 'color': 'gray', 'linewidth': 0.5})
    medianprops: Dict[str, Any] = field(default_factory=lambda:
        {'zorder': 1, 'color': 'gray', 'linewidth': 0.5})

@dataclass(frozen = True)
class LollipopStyleOpts:
    line_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"color":'black', "linewidth":0.5, "linestyle":'--'})
    marker_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"markeredgecolor":"black", "markersize":3, "markeredgewidth":0.4})
    stem_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"color":"black", "linewidth":0.4})
    base_kws: Optional[Dict[str, Any]] = field(default_factory = lambda:
        {"color":"black"})
    marker_rasterized: bool = True
    stem_rasterized: bool = True

@dataclass(frozen = True)
class LoessStyleOpts:
    line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "k", "ls": "--", "lw": 1})
    pval_plot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"alpha": 1, "markersize": 1})
    regr_plot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"alpha": 1, "markersize": 1, "color":"gray"})
    below_color: Optional[str] = "lightgray"
    above_color: Optional[str] = "red"

@dataclass(frozen = True)
class ClusteredHeatmapStyleOpts:
    col_width: float = 0.05
    heatmap_width: float = 0.75
    wspace: float = 0.02
    color_pal: str = "tab20"
    default_color: str = "#EEEEEE"
    palette: Dict[str, Any] = field(default_factory=lambda: None)
    heatmap_kws: Dict[str, Any] = field(default_factory=lambda:
        {"cmap":"bwr", "center":0, "vmin":-10, "vmax":15, "cbar":False, "rasterized":True})
    border_kws: Dict[str, Any] = field(default_factory=lambda:
        {"fill":False, "clip_on":False, "zorder":10, "edgecolor":"black", "linewidth":0.75})
    line_kws: Dict[str, Any] = field(default_factory=lambda:
        {"color":'black', "linewidth":0.75})

@dataclass(frozen = True)
class PWESHeatmapStyleOpts:
    # HEATMAPS #
    heatmap_kws: Dict[str, Any] = field(default_factory=lambda:
        {"square":True, "cmap":'RdBu_r', "vmin":-1, "vmax":1,
         "xticklabels":False, "yticklabels":False, "cbar":False})
    clustered_heatmap_kws: Dict[str, Any] = field(default_factory=lambda:
        {"vmin":-1, "vmax":1, "center":0, "cmap":'RdBu_r',
         "xticklabels":False, "yticklabels":False})
    line_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"color": "k", "ls": "--", "lw": 0.25})

    swarmplot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"size": 2, "linewidth": 0, "edgecolor":None})
    jitterplot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"alpha":0.5, "size":2, "linewidth":0, "edgecolor":None})
    stripplot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"alpha":0.5, "size":2, "linewidth":0, "edgecolor":None})

    boxplot_kws: Optional[Dict[str, Any]] = field(default_factory=lambda:
        {"fliersize":0, "linewidth":0.5, "width":0.5, })
    boxprops: Dict[str, Any] = field(default_factory=lambda:
        {'facecolor': 'none', 'edgecolor': "black"})
    whiskerprops: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'black'})
    capprops: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'black'})
    medianprops: Dict[str, Any] = field(default_factory=lambda:
        {'color': 'black'})

    spine_color: str = "black"
    spine_visible: bool = True
    spine_wspace: float = 0.2
    spine_hspace: float = 0.2
    linewidth: float = 0.5
    edgecolor: str = "black"
    # COLOR BAR #
    colorbar_normalize_kws: Dict[str, Any] = field(default_factory=lambda:
        {"vmin":-1, "vmax":1})
    colorbar_base_kws: Dict[str, Any] = field(default_factory=lambda:
        {"cmap":'RdBu_r'})
    colorbar_linewidth: float = 0.5

@dataclass(frozen = True)
class DomainOpts:
    dom_setting: Dict[str, List[Dict[str, Any]]] = field(default_factory=dict)
    gene: str = field(default_factory=str)
    default_alpha: float = 0.25
    default_color: str = "lightblue"
    gene_col: str = "Gene Symbol"


def resolve_domains(domain, df=None):
    """
    The list of domain spans to shade, from a DomainOpts.

    dom_setting is typed as a mapping of gene name to a list of spans, but was
    being iterated directly, which yields the dict's *keys* -- so any non-empty
    mapping raised TypeError on dom["start"]. This resolves the mapping using
    DomainOpts.gene (previously read by nothing at all), and falls back to the
    only gene present when `gene` is unset and there is exactly one.

    A plain list of spans is also accepted, so callers can pass
    [{'start': 1, 'end': 50}] without wrapping it in a dict.

    Each span is a dict with 'start' and 'end', optionally 'color' and
    'alpha'. Returns an empty list when nothing is configured.
    """
    setting = domain.dom_setting
    if not setting:
        return []
    if isinstance(setting, (list, tuple)):
        return list(setting)
    if not isinstance(setting, dict):
        return []

    if domain.gene and domain.gene in setting:
        return list(setting[domain.gene])
    if len(setting) == 1:
        return list(next(iter(setting.values())))
    if df is not None and domain.gene_col in getattr(df, "columns", []):
        genes = [g for g in df[domain.gene_col].dropna().unique() if g in setting]
        if len(genes) == 1:
            return list(setting[genes[0]])
    return []

@dataclass(frozen = True)
class LegendOpts:
    loc: str = None
    ncol: int = 1
    frameon: bool = False
    mark_size: Optional[int] = 5
    title: Optional[str] = None
    path: Union[str, Path] = "legend"

# Dataclass for outputs
@dataclass(frozen = True)
class OutputOpts:
    figsize: Tuple[float, float] = (4.0, 2.2)
    subplots_kws: Dict[str, Any] = field(default_factory=dict)
    save: bool = False
    path: str = "scatterplot"
    show: bool = True
    dpi: int = 1200
    transparent: bool = True
    rasterize: bool = True
    out_type: str = "svg"

# Defining default options
AXIS = AxisLabelOpts()
SCATTERSTYLE = ScatterStyleOpts()
NEGCONTROL = NegCtrlOpts()
LOLLIPOPSTYLE = LollipopStyleOpts()
BOXPLOTSTYLE = BoxStyleOpts()
JITTERPLOTSTYLE = JitterStyleOpts()
LOESSSTYLE = LoessStyleOpts()
CLUSTEREDHEATMAPSTYLE = ClusteredHeatmapStyleOpts()
PWESHEATMAPSTYLE = PWESHeatmapStyleOpts()
DOMAIN = DomainOpts()
LEGEND = LegendOpts()
OUTPUT = OutputOpts()

def resolve_colors(values, style, palette=None):
    """
    Map category values to colours, for both static and interactive figures.

    ScatterStyleOpts.palette is typed Union[str, List[str]] and defaults to
    the string "colorblind". Zipping values against it, as the code here used
    to, iterates the *characters* of that string, so a hue column produced
    colours 'c', 'o', 'l' and matplotlib rejected them. A string is a seaborn
    palette name; only a list is a list of colours.

    Categories are sorted before assignment so a value keeps its colour
    between a static figure and its interactive twin.
    """
    values = sorted(set(values))
    if palette:
        return {v: palette.get(v) for v in values}
    color_map = getattr(style, 'color_map', None)
    if color_map:
        return {v: color_map.get(v, None) for v in values}

    pal = getattr(style, 'palette', 'colorblind')
    if isinstance(pal, str):
        colors = sns.color_palette(pal, n_colors=len(values))
    else:
        colors = list(pal)
    if not colors:
        colors = sns.color_palette('colorblind', n_colors=max(len(values), 1))
    hexes = [c if isinstance(c, str) else mcolors.to_hex(c) for c in colors]
    return {v: hexes[i % len(hexes)] for i, v in enumerate(values)}


def load_frame(df):
    """
    Accept a DataFrame or a path to a .csv, and return a copy.

    Only two figure_plot functions did this; the rest required a DataFrame and
    several mutated the caller's frame in place.
    """
    if isinstance(df, (str, Path)):
        return pd.read_csv(df)
    return df.copy()


def apply_axis_limits(ax, axis):
    """
    Apply whichever of xlim/ylim/xticks/yticks are set.

    Several functions called ax.set(ylim=..., yticks=...) unconditionally.
    Matplotlib's set_ticks(None) raises TypeError, so those functions crashed
    on the default AxisLabelOpts, where yticks is None.
    """
    if axis.xlim is not None:
        ax.set_xlim(*axis.xlim)
    if axis.ylim is not None:
        ax.set_ylim(*axis.ylim)
    if axis.xticks is not None:
        ax.set_xticks(axis.xticks)
    if axis.yticks is not None:
        ax.set_yticks(axis.yticks)


def figure_path(output, suffix='', ext=None):
    """
    Where a figure should be written, with its parent directory created.

    suffix distinguishes the several files a single function may emit, such as
    a plot and its separate legend.
    """
    ext = output.out_type if ext is None else ext
    path = Path(f"{str(output.path)}{suffix}.{ext}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def save_figure(fig, output, suffix='', **savefig_kws):
    """
    Write a figure out, but only when output.save is set.

    boxplot_figure, jitterbox_kdeplot_figure and kmeans_clustered_heatmap all
    saved unconditionally, so output.save did nothing for them and simply
    drawing a figure wrote files to disk.
    """
    if not output.save:
        return None
    path = figure_path(output, suffix)
    kws = {'dpi': output.dpi, 'format': output.out_type,
           'transparent': output.transparent}
    kws.update(savefig_kws)
    fig.savefig(path, **kws)
    return path


def finish_figure(fig, output, suffix='', **savefig_kws):
    """Save if asked, show if asked, and close otherwise. Returns the figure."""
    save_figure(fig, output, suffix=suffix, **savefig_kws)
    if output.show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def _load_figure_font(size=6, filename='Arial.ttf'):
    """
    Font used for every figure_plot label.

    This was FontProperties(fname='Arial.ttf'), a path resolved against the
    caller's working directory. Matplotlib does not open the file until it
    renders, so a figure only failed at savefig/tight_layout time, and only
    for users who happened not to have Arial.ttf in the directory they ran
    from. The lookup order below keeps that file authoritative when present
    and degrades to a system font instead of raising.
    """
    for candidate in (Path.cwd() / filename, Path(__file__).parent / filename):
        if candidate.is_file():
            return fm.FontProperties(fname=str(candidate), size=size)
    try:
        found = fm.findfont(fm.FontProperties(family='Arial'),
                            fallback_to_default=False)
        return fm.FontProperties(fname=found, size=size)
    except Exception:
        return fm.FontProperties(family='sans-serif', size=size)


arial_font6 = _load_figure_font(size=6)
 
