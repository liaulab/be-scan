# Publication Figure Functions: `figure_plot`

`be_scan.figure_plot` produces figures intended for a paper. It differs from
{doc}`plot` in two ways.

**Options arrive as frozen dataclasses** rather than loose keyword arguments,
so one style can be defined once and reused across every panel of a figure:

```python
from dataclasses import replace
from be_scan.figure_plot import AxisLabelOpts, OutputOpts, ScatterStyleOpts

AXIS  = AxisLabelOpts(xlabel="residue", ylabel="log2 fold change")
STYLE = ScatterStyleOpts(hue_col="muttype",
                         color_map={"Missense": "#fb8072", "Silent": "#d9d9d9"})
OUT   = OutputOpts(figsize=(5.0, 2.2), out_type="svg", dpi=600, save=True)

TALL  = replace(OUT, figsize=(5.0, 3.2))   # a variant; the originals are frozen
```

**Every figure function has an interactive counterpart.** The twin takes the
same arguments, is styled from the same option objects, and writes its `.html`
beside the static file:

```python
scatterplot_figure(df, ["score"], "pos", output=OUT)              # figures/screen.svg
interactive_scatterplot_figure(df, ["score"], "pos", output=OUT)  # figures/screen.html
```

| static | interactive |
| --- | --- |
| `scatterplot_figure` | `interactive_scatterplot_figure` |
| `correlation_scatterplot_figure` | `interactive_correlation_scatterplot_figure` |
| `boxplot_figure` | `interactive_boxplot_figure` |
| `jitterbox_kdeplot_figure` | `interactive_jitterbox_kdeplot_figure` |
| `lollipop_figure` | `interactive_lollipop_figure` |
| `lollipop_split_figure` | `interactive_lollipop_split_figure` |
| `domain_plot_figure` | `interactive_domain_plot_figure` |
| `kmeans_clustered_heatmap` | `interactive_kmeans_clustered_heatmap` |
| `loess_figure` | `interactive_loess_figure` |
| `pwes_histogram` | `interactive_pwes_histogram` |
| `pwes_boxplots` | `interactive_pwes_boxplots` |
| `pwes_scatterplots` | `interactive_pwes_scatterplots` |
| `pwes_heatmap` | `interactive_pwes_heatmap` |
| `pwes_clustermap` | `interactive_pwes_clustermap` |
| `pwes_heatmap_colorbar` | `interactive_pwes_heatmap_colorbar` |

## How the two are kept consistent

The plotly layout is derived from the live matplotlib rcParams the static
figures use, and marker and line styling is read from the same dataclass
fields, so restyling one restyles the other. Three conversions make the
renderings actually match:

- matplotlib's scatter `s` is an **area** in points squared while plotly's
  marker size is a **diameter**, so the conversion is a square root
- the style dataclasses use matplotlib colour shorthands such as `"k"`, which
  plotly rejects, so every colour is normalised first
- `OutputOpts.figsize` is in inches for print, typically 4 × 2.2 at 6 pt,
  which is illegible in a browser. Each interactive function takes a `scale`
  that multiplies pixel size and font size by the same factor, preserving
  proportions exactly. `scale=1.0` reproduces the print geometry.

## File size

`plotly.js` is embedded by default, so each `.html` is self-contained
(≈3.7 MB) and opens with no network connection. For many figures at once:

```python
from be_scan.figure_plot import figure_interactive_style as ivs
ivs.INCLUDE_PLOTLYJS = "cdn"     # ≈40 KB per file, requires a connection
```

## Domain maps

`domain_plot_figure` draws a protein's domains as a bar where every residue
has the same width. The axes fill the figure, so at the default 4 in width the
map lines up under a scatterplot or lollipop of the same width.

```python
gene_domains = [
    {"description": "IDR1", "start": 1,   "end": 175, "color": "#fb8072"},
    {"description": "QLQ",  "start": 176, "end": 207, "color": "#80b1d3"},
    # ...
]
domain_plot_figure(gene_domains, gene_len=1590, label_angle=45,
                   output=replace(DOMAINPLOTOUTPUT, save=True, path="figures/smarca4"))
```

A domain is filled with its `color` only when that is a hex code, otherwise
with `DomainPlotStyleOpts.default_color` (`#b7b7b7`). Boxes have a 0.5 pt black
outline and rounded corners over a 0.5 pt black baseline. Labels are 6 pt
Arial, sit under their domain, and are pushed sideways just enough that none
overlap; `label_angle` rotates them. `DOMAINPLOTOUTPUT` is the scatterplot
width at a tenth of its height, and the file is saved with a tight bounding
box so the labels hanging below are kept.

## Option objects

| dataclass | controls |
| --- | --- |
| `AxisLabelOpts` | title, labels, limits, ticks, x window |
| `ScatterStyleOpts` | palette, hue and marker columns, per-layer keywords, transparency threshold |
| `BoxStyleOpts`, `JitterStyleOpts`, `LollipopStyleOpts`, `LoessStyleOpts`, `ClusteredHeatmapStyleOpts`, `PWESHeatmapStyleOpts` | per-plot styling |
| `DomainOpts` | shaded protein domains; accepts a list of spans or a mapping of gene to spans |
| `DomainPlotStyleOpts` | domain map fill, outline, baseline, box height, corner rounding, label spacing |
| `LegendOpts` | placement, columns, and an optional separate legend file |
| `NegCtrlOpts` | negative-control normalisation and ± SD guide lines |
| `OutputOpts` | figure size, dpi, format, transparency, whether to save and show |

See {doc}`notebooks/04_figure_plot` for a worked example of each.

```{eval-rst}
.. automodule:: be_scan.figure_plot
    :members:
    :undoc-members:
    :imported-members:
```
