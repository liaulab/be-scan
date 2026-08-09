# Changelog

## 3.0.0

The guiding constraint for this release was that existing code should keep
working. Every new parameter is appended to its signature and defaults to the
previous behavior, and a test suite compares all 40 public functions against a
snapshot of the 2.0.0 API on every run, failing if a parameter is renamed,
reordered, removed, or given a different default.

Where behavior did change, it is listed under *Behavior changes* below.

### New

**Replicate significance.** For every condition with more than one replicate,
`log_transform` now adds six columns alongside the existing `_avg` and
`_stdev`:

| column | meaning |
| --- | --- |
| `_pval` | two-sided one-sample t-test of the replicate log-fold-changes against 0 |
| `_neglogp` | −log10 of it |
| `_pval_mod` | the same test with variance shrunk toward a prior fitted across the library |
| `_neglogp_mod` | −log10 of it |
| `_qval` | Benjamini-Hochberg FDR of `_pval` |
| `_qval_mod` | Benjamini-Hochberg FDR of `_pval_mod` |

Single-replicate conditions get none of these, since a p-value across one
measurement is not defined. No new parameters are required.

**With three replicates, prefer the moderated columns.** The ordinary t-test
ranks largely on how lucky each guide's variance estimate happened to be. On
the FOXA1 example screen its top hit has a log-fold-change of 0.155 with a
standard deviation of 0.003, and a non-targeting control ranks fifth of 1291.
After FDR correction the ordinary test finds nothing significant anywhere in
that screen; the moderated test finds 109 guides, led by the dependency
controls. Read either alongside `_avg` — a small p-value on a negligible
effect is still a negligible effect.

Guides absent from every replicate have identical log-fold-changes and no
variance. `scipy` reports p=0 for that, which is −log10(p) ≈ 308 and would
top every volcano plot, so those rows are `NaN`.

**Tunable spacer length.** `guide_len` is added to `generate_library`,
`annotate`, `count_reads`, `batch_process` and `GeneForCRISPR`, defaulting to
20. The searched sequence is `guide_len + len(PAM)` rather than a fixed 23, so
non-3 nt PAMs work too.

**Interactive figures.** Every function in `figure_plot` has a plotly
counterpart taking the same arguments and styled from the same option objects,
writing its `.html` beside the static file. See {doc}`figure_plot`.

**Consistent plotting options.** Every function in `be_scan.plot` now accepts
the same set: negative-control normalisation, value and category filters, hue,
subset highlighting, labels and title, output controls including `dpi`. Each
returns `(fig, axes)`.

### Behavior changes

These alter output. All are corrections.

- **`window` is genuinely tunable.** It was documented as tunable but only
  correct at its default; of every window on a 20 nt guide, one produced
  answers matching an independent reference implementation, and it was `[4,8]`.
  The rest either raised or silently shifted residue numbering. Any window with
  `start >= 3` and `end <= 18` is now correct. Output at `[4,8]` is unchanged
  apart from the two `windowend_pos` corrections below.

- **`windowend_pos`** was decided by the base *after* the editing window rather
  than its last base. On the FOXA1 example this changed 2 guides of 992: one
  all-exonic window that reported `-1`, and one all-intronic window that
  reported a position. One downstream `_pos` value changed as a result. Across
  8,946 guide × editor comparisons against 2.0.0, nothing else differs.

- **`scatterplot`'s filters now filter.** `filter_val` and `filter_params` had
  no effect on the plotted points: the plot call captured a reference to the
  pre-filter frame while the annotation step used the filtered one, so the two
  disagreed. A 300-row frame filtered to `["Missense"]` plotted all 300 rows
  instead of 63.

- **`corr_jointplot`'s legend** was attached to the y-marginal histogram rather
  than the joint axes, producing an empty legend box beside the marginal. It is
  now drawn on the joint axes, and only when `include_hue` is set.

- **`pwes_clustering`'s `out_dir`** is joined as a path rather than
  concatenated. `out_dir='results'` previously created an empty `results/` and
  wrote every file to `resultstest_*.csv` beside it. A trailing separator, the
  documented usage, is unaffected.

- **`output.save` is honoured throughout `figure_plot`.** `boxplot_figure`,
  `jitterbox_kdeplot_figure`, `kmeans_clustered_heatmap`, five of the PWES
  plots, and `scatterplot_figure`'s separate legend file all saved
  unconditionally, so drawing a figure wrote files to disk.

- **`kmeans_clustered_heatmap`** no longer calls `plt.show()` regardless of
  `output.show`, and returns its figure.

### Fixed

Each of these was reachable on a default or documented path.

- `calc_controls` tested `isinstance(df_conds, ...)` before `df_conds` was
  assigned, so every call raised `UnboundLocalError` and the statistics branch
  of `batch_process` was dead code.
- `correlation_scatterplot_figure` closed its legend figure outside the block
  that creates it, raising `NameError` whenever `output.save` was `False`.
- `scatterplot_figure` zipped hue values against `style.palette`, which
  defaults to the string `"colorblind"`, iterating that string's characters and
  handing matplotlib the colours `'c'`, `'o'`, `'l'`. Any call with a hue
  column and the default palette raised `ValueError`.
- `figure_plot` built its font from the relative path `Arial.ttf`, resolved at
  render time against the caller's working directory, so every figure raised
  `FileNotFoundError` unless run from a directory containing that file.
- `DomainOpts.dom_setting` is a mapping of gene to spans but was iterated
  directly, yielding dict keys, so any non-empty setting raised `TypeError` in
  the scatterplot, lollipop, loess and PWES figures.
- `scatterplot_figure` required a `Gene Symbol` column on every call, even with
  no domains configured.
- `plot.scatterplot` called IPython's `display()` unconditionally, so
  `annot_top` and `annot_abs` raised `NameError` outside a notebook; `annot_abs`
  was also missing from the validation gate, making it unreachable.
- `boxplot_figure` and `jitterbox_kdeplot_figure` called
  `ax.set(yticks=None)`, which raises, so neither worked with the default
  `AxisLabelOpts`.
- `pwes_clustering` raised `TypeError` with its default `out_dir=None`, and
  mutated the caller's frame.
- `pwes_histogram` and `pwes_scatterplots` read style fields that exist only on
  `PWESHeatmapStyleOpts` while defaulting to `ScatterStyleOpts`, raising
  `AttributeError`; `pwes_heatmap`'s `blackout` assumed string labels and raised
  on the single-subunit path; PWES colour lists now cycle rather than raising
  `IndexError` when shorter than the cluster count.
- `fill_in_dna_exons` indexed `exons[exon_i - 1]` into the last exon and
  `exons[exon_i + 1]` past the end without bounds checks.
- `annotate` ignored its own `seq_col`, `window_start_col` and `window_end_col`
  parameters at four call sites, so renamed columns raised `KeyError` or were
  silently the wrong column. `gene_col` is added for the same reason.
- `count_reads` named `annotated_lib` in an error validating the sample sheet.
- Plot functions no longer mutate the caller's DataFrame during
  negative-control normalisation.

### Testing

The suite covers every public function. The annotation path is cross-checked
against a second implementation sharing none of its arithmetic: `be_scan`
derives mutations through per-guide frame tracking and codon-window slicing,
while the reference works in absolute CDS coordinates and retranslates the
whole sequence, where a codon split across a junction was never split. The two
agree on 33,801 mutation calls across ABE, CBE and dual editors on real FOXA1
input, and on synthetic genes built with junctions at all three codon phases.

## 2.0.0

Baseline release, tagged to mark the behavior existing code depends on. No
functional change from the previously published version.
