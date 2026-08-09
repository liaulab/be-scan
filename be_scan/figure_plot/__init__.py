
from be_scan.figure_plot.figure_classes import *

from be_scan.figure_plot.figure_scatterplot import scatterplot_figure
from be_scan.figure_plot.figure_correlation_scatterplot import correlation_scatterplot_figure
from be_scan.figure_plot.figure_boxplot import boxplot_figure
from be_scan.figure_plot.figure_jitterplot import jitterbox_kdeplot_figure
from be_scan.figure_plot.figure_lollipop import lollipop_figure, lollipop_split_figure

from be_scan.figure_plot.figure_clustered_heatmap_main import kmeans_clustered_heatmap
from be_scan.figure_plot.figure_clustered_heatmap_helpers import kmeans_elbow, kmeans_pca_scatterplot, hex_to_rgb

from be_scan.figure_plot.figure_loess_main import loess_smoothing
from be_scan.figure_plot.figure_loess_helpers import interp_methods
from be_scan.figure_plot.figure_loess_plot import loess_figure

from be_scan.figure_plot.figure_pwes_main import pwes_clustering
from be_scan.figure_plot.figure_pwes_plot import pwes_histogram, pwes_boxplots, pwes_scatterplots, pwes_heatmap, pwes_clustermap, pwes_heatmap_colorbar
from be_scan.figure_plot.figure_pwes_helpers import fillgaps

# INTERACTIVE COUNTERPARTS #
# Each mirrors its static twin's signature and is driven by the same option
# objects, and writes its .html beside the static file.
from be_scan.figure_plot.figure_interactive import (
    interactive_scatterplot_figure,
    interactive_correlation_scatterplot_figure,
    interactive_boxplot_figure,
    interactive_jitterbox_kdeplot_figure,
    interactive_lollipop_figure,
    interactive_lollipop_split_figure,
    interactive_kmeans_clustered_heatmap,
    interactive_loess_figure,
    interactive_pwes_histogram,
    interactive_pwes_boxplots,
    interactive_pwes_scatterplots,
    interactive_pwes_heatmap,
    interactive_pwes_clustermap,
    interactive_pwes_heatmap_colorbar,
)

# figure_classes re-exports pandas, numpy, matplotlib, seaborn, scipy and
# sklearn names for the convenience of the modules that do `import *` from it.
# Without an explicit __all__ those leak into be_scan.figure_plot, and autodoc
# then pulls third-party docstrings into the API reference.
__all__ = [
    # option dataclasses and their defaults
    "AxisLabelOpts", "NegCtrlOpts", "ScatterStyleOpts", "BoxStyleOpts",
    "JitterStyleOpts", "LollipopStyleOpts", "LoessStyleOpts",
    "ClusteredHeatmapStyleOpts", "PWESHeatmapStyleOpts", "DomainOpts",
    "LegendOpts", "OutputOpts",
    "AXIS", "NEGCONTROL", "SCATTERSTYLE", "BOXPLOTSTYLE", "JITTERPLOTSTYLE",
    "LOLLIPOPSTYLE", "LOESSSTYLE", "CLUSTEREDHEATMAPSTYLE", "PWESHEATMAPSTYLE",
    "DOMAIN", "LEGEND", "OUTPUT",
    # shared helpers
    "resolve_colors", "resolve_domains", "load_frame", "figure_path",
    "save_figure", "finish_figure", "apply_axis_limits",
    # static figures
    "scatterplot_figure", "correlation_scatterplot_figure", "boxplot_figure",
    "jitterbox_kdeplot_figure", "lollipop_figure", "lollipop_split_figure",
    "kmeans_clustered_heatmap", "kmeans_elbow", "kmeans_pca_scatterplot",
    "hex_to_rgb", "loess_smoothing", "interp_methods", "loess_figure",
    "pwes_clustering", "pwes_histogram", "pwes_boxplots", "pwes_scatterplots",
    "pwes_heatmap", "pwes_clustermap", "pwes_heatmap_colorbar", "fillgaps",
    # interactive counterparts
    "interactive_scatterplot_figure",
    "interactive_correlation_scatterplot_figure",
    "interactive_boxplot_figure", "interactive_jitterbox_kdeplot_figure",
    "interactive_lollipop_figure", "interactive_lollipop_split_figure",
    "interactive_kmeans_clustered_heatmap", "interactive_loess_figure",
    "interactive_pwes_histogram", "interactive_pwes_boxplots",
    "interactive_pwes_scatterplots", "interactive_pwes_heatmap",
    "interactive_pwes_clustermap", "interactive_pwes_heatmap_colorbar",
]
