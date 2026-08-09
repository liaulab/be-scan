# BE-SCAN

**B**ase-**E**diting **S**creening **C**omputational **AN**alysis — a Python
package for designing, analysing and plotting CRISPR base editor and
mutagenesis screens.

Authors: Calvin XY Hu, Simon P Shen, Christine Lee, Brian B Liau
· Liau Lab, Harvard University

Documentation: <https://be-scan.readthedocs.io>

## Install

```console
pip install git+https://github.com/liaulab/be-scan.git
```

Takes about five minutes. Python 3.7 or newer.

## The four modules

| module | what it does | |
| --- | --- | --- |
| `be_scan.sgrna` | design a library from a gene sequence and annotate every guide with the residues it could change | [readme](be_scan/sgrna/README.md) |
| `be_scan.analysis` | count reads from FASTQ, normalise, compare conditions, test replicates | [readme](be_scan/analysis/README.md) |
| `be_scan.plot` | exploratory plots, static and interactive | [readme](be_scan/plot/README.md) |
| `be_scan.figure_plot` | publication figures, each with an interactive counterpart | [readme](be_scan/figure_plot/README.md) |

## A minimal run

```python
from be_scan.sgrna import design_library
from be_scan.analysis import batch_process
from be_scan.plot import scatterplot

# 1. design and annotate a library from a gene .fasta
design_library(gene_filepath="FOXA1.fasta", gene_name="FOXA1",
               cas_type="SpG", output_name="library.csv")

# 2. count reads and score guides
batch_process(sample_sheet="samples.csv", annotated_lib="library.csv",
              comparisons="comparisons.csv", controls=["t0"])

# 3. look at the result
scatterplot(df_filepath="conditions.csv", comparisons=["sorted-unsorted"],
            x_column="AtoG_pos")
```

## Tutorials

Four notebooks under `docs/source/notebooks/`, each running end to end on the
included FOXA1 example data:

1. `01_design` — from a gene `.fasta` to an annotated library
2. `02_analysis` — from FASTQ to per-guide scores
3. `03_plot` — exploring the result
4. `04_figure_plot` — figures for a paper

## Version 3.0.0

Existing code keeps working. Every new parameter is appended to its signature
and defaults to the previous behavior, and the test suite compares all 40
public functions against a snapshot of the 2.0.0 API on every run, failing if
a parameter is renamed, reordered, removed, or given a different default.

The handful of places where output genuinely changed are listed in
[CHANGELOG.md](CHANGELOG.md); all are corrections.

Headline changes:

- **replicate significance** — `log_transform` reports p-values, variance-shrunk
  p-values and FDR-corrected q-values, with no new parameters
- **`window` is genuinely tunable** — it was documented as such but only correct
  at its default
- **`guide_len` is a parameter** — 20 nt spacers were assumed in six places
- **interactive figures** — every `figure_plot` function has a plotly twin
  writing its `.html` beside the static file
- **consistent plotting options** — every `be_scan.plot` function takes the same
  filters, colours, highlighting and output controls

## Development

```console
git clone https://github.com/liaulab/be-scan.git
cd be_scan
pip install -e ".[test,docs]"
pytest
```

The suite is ~700 tests. The annotation path is cross-checked against a second
implementation sharing none of its arithmetic; see `tests/sgrna_test/`.

### Building the documentation

```console
pip install -e ".[docs]"
sphinx-autobuild docs/source docs/build/html/
```

The tutorial notebooks are committed with their outputs, and `conf.py` sets
`nbsphinx_execute = 'never'`, so a docs build needs neither a kernel nor the
example data. Rendering the notebooks' markdown cells does need **pandoc**,
which is a system package rather than a Python one:

```console
conda install -c conda-forge pandoc     # or: brew install pandoc
```

### Read the Docs

If `/en/latest/` is stale while `/en/main/` is current, the cause is that
`latest` is not "the most recent build" — it is a permanent alias for whatever
Read the Docs has recorded as the project's **default branch**. A project set
up when the branch was `master` keeps pointing there after a rename, so
`latest` rebuilds a branch that no longer receives commits, while the
separately-created `main` version tracks the real one.

To fix it, in the Read the Docs dashboard:

1. **Admin → Settings → Default branch** → set to `main`, and save.
2. **Admin → Versions** → confirm `latest` is Active.
3. **Builds → Build version: `latest`** → trigger one build. `/en/latest/` now
   matches `/en/main/`.
4. **Admin → Settings → Default version** → set to `latest`, so the bare
   `be-scan.readthedocs.io` lands there.
5. **Admin → Redirects** → add `/en/main/` → `/en/latest/`, so the URL already
   shared elsewhere keeps working, then deactivate the `main` version rather
   than maintaining two identical copies.
6. **Admin → Integrations** → check the GitHub webhook exists and shows recent
   deliveries. Without it, pushes never trigger a build regardless of the
   branch setting.

Once tags are pushed, activating them under **Admin → Versions** gives a
`stable` alias pinned to the highest tag, so 2.0.0 and 3.0.0 docs stay
browsable side by side.

## Citation

If BE-SCAN is useful in your work, please cite the Liau Lab.

## License

MIT — see [LICENSE.txt](LICENSE.txt).
