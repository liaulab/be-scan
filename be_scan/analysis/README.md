# `be_scan.analysis` — FASTQ to scores

Counts reads, normalises them, compares conditions, and tests replicates for
significance.

## Public functions

| function | purpose |
| --- | --- |
| `count_reads` | extract spacers from FASTQ and match them to the library |
| `log_transform` | reads per million, log2, subtract control, summarise replicates |
| `compare_conds` | subtract one condition group's average from another's |
| `calc_controls` | mean and spread of the negative controls |
| `batch_process` | run all four in order |

## The sample sheet

Drives everything. Required columns:

| column | meaning |
| --- | --- |
| `fastq_file` | reads for this sample |
| `condition` | this sample's own name, unique per row |
| `agg_conditions` | which group it replicates |
| `counts_file`, `noncounts_file`, `stats_file` | per-sample outputs |

`condition` versus `agg_conditions` is the distinction that matters:
`sorted-top-1`, `-2` and `-3` are three conditions belonging to one
`sorted_top` group, and that grouping is what replicate statistics are
computed across.

## Quick start

```python
from be_scan.analysis import batch_process

batch_process(
    sample_sheet  = "samples.csv",
    annotated_lib = "library.csv",
    comparisons   = "comparisons.csv",
    controls      = ["t0"],
    guide_len     = 20,
)
```

## Read counting

`count_reads` searches each read for `KEY`, takes the spacer that follows, and
stops at `KEY_REV`.

| parameter | default | meaning |
| --- | --- | --- |
| `KEY` | `CGAAACACC` | end of the hU6 promoter |
| `KEY_REV` | `GTTTGAGA` | start of the scaffold; empty to take a fixed length |
| `KEY_INTERVAL` | `(10, 80)` | where in the read to look for `KEY` |
| `dont_trim_G` | `False` | keep a leading G on spacers one base too long |
| `guide_len` | `20` | expected spacer length |

## Significance

For every condition with **more than one replicate**, `log_transform` adds:

| column | meaning |
| --- | --- |
| `_pval` | two-sided one-sample t-test of the replicate LFCs against 0 |
| `_neglogp` | −log10 of it |
| `_pval_mod` | the same test with variance shrunk toward a prior fitted across the library |
| `_neglogp_mod` | −log10 of it |
| `_qval` | Benjamini-Hochberg FDR of `_pval` |
| `_qval_mod` | Benjamini-Hochberg FDR of `_pval_mod` |

Single-replicate conditions get none, since a p-value across one measurement
is not defined.

### Which to use

**With two or three replicates, prefer the moderated columns.** The ordinary
t-test ranks largely on how lucky each guide's variance estimate happened to
be. On the FOXA1 example screen its top hit has a log-fold-change of 0.155
with a standard deviation of 0.003, and a non-targeting control ranks fifth of
1291. After FDR correction the ordinary test finds nothing significant
anywhere in that screen; the moderated test finds 109, led by the dependency
controls.

Read either alongside `_avg`. A small p-value on a negligible effect is still
a negligible effect.

### Guides with no variance

An sgRNA absent from every replicate has identical log-fold-changes and no
spread. `scipy` reports p=0 for that, which is −log10(p) ≈ 308 and would sit
at the top of every volcano plot. Those rows are `NaN` instead.

## Normalisation

`controls=["t0"]` averages the named samples and subtracts them, producing
`_LFCminusControl` columns. `controls=[]` skips that step and the suffix is
`_LFC`. Whichever you choose, use it consistently — `compare_conds` and
`calc_controls` pick their columns from the same `controls` argument.
