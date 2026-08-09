# `be_scan.plot` — exploratory plots

The quick-look layer: ordinary keyword arguments, sensible defaults,
matplotlib and plotly side by side. Use it while working out what a screen
says; use {doc}`../figure_plot` for figures that go into a paper.

## Public functions

| static | interactive |
| --- | --- |
| `scatterplot` | `interactive_scatter` |
| `boxplot` | `interactive_boxplot` |
| `corr_jointplot` | `interactive_corr_jointplot` |
| `corr_heatmap` | `interactive_corr_heatmap` |
| `sgrna_enrichment` | — |
| `loess_smoothing` (+ `plot_loess`) | — |

## The shared options

Every function accepts the same set. All default to off, so adding one only
changes the plot when you ask.

| group | parameters |
| --- | --- |
| data | `df_filepath` — a DataFrame **or** a path to a csv |
| negative control | `neg_ctrl`, `neg_ctrl_col`, `neg_ctrl_conditions` |
| filter by value | `filter_val`, `val_cols`, `val_min` |
| filter by category | `filter_params`, `params_cols`, `params_conditions` |
| colour | `include_hue`, `hue_col`, `pal` |
| highlight | `highlight`, `highlight_col`, `highlight_vals` |
| annotate | `annot`, `annot_label`, `annot_top`, `annot_cutoff`, `annot_abs` |
| labels | `xlab`, `ylab`, `title` |
| output | `savefig`, `show`, `out_name`, `out_type`, `out_dir`, `dpi` |

Each returns `(fig, axes)`, so you can keep customising after the call.

Before 3.0.0 these were scattered: `scatterplot` could filter, colour and
annotate while `corr_heatmap` could do none of it, and only
`sgrna_enrichment` could highlight a subset.

## Quick start

```python
import be_scan.plot as bsp

fig, axes = bsp.scatterplot(
    df_filepath = "conditions.csv",
    comparisons = ["sorted-unsorted"],
    x_column    = "AtoG_pos",

    include_hue = True, hue_col = "AtoG_muttype",
    filter_params = True, params_cols = ["AtoG_muttype"],
    params_conditions = [["Missense", "Nonsense", "Silent"]],
    neg_ctrl = True, neg_ctrl_col = "gene", neg_ctrl_conditions = ["Control"],
    domains  = [{"start": 168, "end": 258}],

    xlab = "residue", ylab = "log2 fold change",
)
```

## Annotation modes

| parameter | labels |
| --- | --- |
| `annot_top=n` | the n highest |
| `annot_top=-n` | the n lowest |
| `annot_abs=n` | the n largest by absolute value |
| `annot_cutoff=v` | everything beyond ±v |

## Negative controls

`neg_ctrl=True` subtracts the mean of the control guides from every
comparison, and draws ±`n_sd` × their standard deviation as dashed lines.
`n_sd` defaults to 2. The caller's DataFrame is never modified.

## Two defaults that look wrong and are kept anyway

Correcting either would silently change results for existing code, so both are
left as they are:

- `corr_jointplot`'s `hue_col` defaults to `'CtoT_muttype'` where
  `scatterplot` uses `''`
- `sgrna_enrichment`'s `out_name` defaults to `'scatterplot'`

Pass them explicitly.

## LOESS

`loess_smoothing` fits a curve along the protein and tests it against
randomised data, so it finds *regions* behaving unusually rather than
individual guides. It returns a results dict; pass it to `plot_loess`.

```python
from be_scan.plot.loess_smoothing import loess_smoothing, plot_loess

results = loess_smoothing(df=coding, x_column="AtoG_pos",
                          comparisons=["sorted-unsorted"],
                          span=20, n_repeats=1000)
plot_loess(results, gene="FOXA1", domains=[{"start": 168, "end": 258}])
```

`n_repeats` drives the runtime: 1,000 takes seconds, 10,000 a couple of
minutes.
