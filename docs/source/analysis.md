# Library and Screen Analysis Functions: `analysis`

```{eval-rst}
.. automodule:: be_scan.analysis
    :members:
    :undoc-members:
    :imported-members:
```

## count_reads()

This function counts the reads in a FASTQ file and assign them to a reference sgRNA set {py:func}`be_scan.analysis.count_reads`

```{eval-rst}
.. autofunction:: be_scan.analysis.count_reads
```

## merge_and_norm()

This function aggregate raw counts and perform log2-transform and t0 normalization {py:func}`be_scan.analysis.merge_and_norm`

```{eval-rst}
.. autofunction:: be_scan.analysis.merge_and_norm
```

## average_reps()

This function for a set of conditions, average the replicates and export it to csv {py:func}`be_scan.analysis.average_reps`

```{eval-rst}
.. autofunction:: be_scan.analysis.average_reps
