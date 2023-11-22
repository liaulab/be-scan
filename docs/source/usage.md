# Usage

## Installation

To install be_scan: (this doesn't work yet)

```console
(.venv) $ pip install be_scan
```

## Genomic Functions

This function converts a DNA sequence to and amino acid sequence {py:func}`be_scan.sgrna.DNA_to_AA`

```{eval-rst}
.. function:: be_scan.sgrna.DNA_to_AA(seq)

    can write in docs manually

    :seq kind: string of base pairs (ACGT)
```

This function takes the reverse complement of a DNA sequence {py:func}`be_scan.sgrna.rev_complement`

```{eval-rst}
.. autofunction:: be_scan.sgrna.rev_complement

    hi

    this requires the docs under the function to be in a particular format tho
```

This function takes in a dataframe from count_reads, performs normalization, and then plots the data for each condition to reveal which guides are enriched {py:func}`be_scan.plot.plot_scatterplot`

```{eval-rst}
.. autofunction:: be_scan.plot.plot_scatterplot
```

This function takes in a dataframe from count_reads, and plots a heatmap showing correlation between all listed conditions {py:func}`be_scan.plot.plot_corr_scatterplot`

```{eval-rst}
.. autofunction:: be_scan.plot.plot_corr_scatterplot
```

This function takes in a dataframe from count_reads, and plots a heatmap showing correlation between all listed conditions {py:func}`be_scan.plot.plot_corr_heatmap`

```{eval-rst}
.. autofunction:: be_scan.plot.plot_corr_heatmap
```
