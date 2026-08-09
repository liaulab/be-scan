# be_scan

**B**ase-**E**diting **S**creening **C**omputational **AN**alysis — a Python
package for designing, analysing and plotting CRISPR base editor and
mutagenesis screens.

The package is four modules, following the shape of a screen:

| module | what it does |
| --- | --- |
| {doc}`sgrna` | design a library from a gene sequence, and annotate every guide with the residues it could change |
| {doc}`analysis` | count reads from FASTQ, normalise, compare conditions, and test replicates for significance |
| {doc}`plot` | exploratory plots, static and interactive |
| {doc}`figure_plot` | publication figures, each with an interactive counterpart |

## Getting started

New to the package: work through the tutorials in order. Each runs end to end
on the included FOXA1 example data.

1. {doc}`notebooks/01_design` — from a gene `.fasta` to an annotated library
2. {doc}`notebooks/02_analysis` — from FASTQ to per-guide scores
3. {doc}`notebooks/03_plot` — exploring the result
4. {doc}`notebooks/04_figure_plot` — figures for a paper

Installation instructions are in {doc}`install`.

## A minimal run

```python
from be_scan.sgrna import design_library
from be_scan.analysis import batch_process
from be_scan.plot import scatterplot

# 1. design and annotate a library
design_library(gene_filepath="FOXA1.fasta", gene_name="FOXA1",
               cas_type="SpG", output_name="library.csv")

# 2. count reads and score guides
batch_process(sample_sheet="samples.csv", annotated_lib="library.csv",
              comparisons="comparisons.csv", controls=["t0"])

# 3. look at the result
scatterplot(df_filepath="conditions.csv", comparisons=["sorted-unsorted"],
            x_column="AtoG_pos")
```

```{toctree}
:caption: 'Getting started'
:maxdepth: 2

install
notebooks/01_design
notebooks/02_analysis
notebooks/03_plot
notebooks/04_figure_plot
```

```{toctree}
:caption: 'API reference'
:maxdepth: 2

sgrna
analysis
plot
figure_plot
```

```{toctree}
:caption: 'Project'
:maxdepth: 1

changelog
```

# Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
