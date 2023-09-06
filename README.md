# BIGSCAM (tentatively)
Base-editing Indels and Genomic Screen with CRISPR And More \
BioInformatics, Genomic Screening, and Computational Analysis for Mutational data \
Authors: SPS, XYH \

Description: We are building a package for all of the lab's computational needs. ^_^ This package is compatible with CRISPOR generated data. 


# Documentation


## sgRNA

### bigscam.sgrna.find_all()

Generates a list of all single guide RNAs that can be edited by this gene
```
bigscam.sgrna.find_all_BE(gene, PAM, edit, window)
bigscam.sgrna.find_all_DSB(gene, PAM, residue)
```
Parameters: 
- gene: a fasta file generated from
- PAM: a recognition motif specific the the Cas9 type
- edit: base edit (ABE A to T, CBE C to G, )
- residue: residue where Cas9 cuts
- window: window of edit inclusive (4 to 8, 4 to 7)

```
bigscam.sgrna.find_all_CRISPOR(gene, species, type)
```
Queries CRISPOR for a list of guides
Parameters: 
- gene: name or sequence
- species: 
- type: Cas9 type

### bigscam.sgrna.annotate()
```
bigscam.sgrna.annotate_BE(sgrnas, gene)
bigscam.sgrna.annotate_DSB(sgrnas, gene)
bigscam.sgrna.annotate_CRISPOR(sgrnas, gene)
```
Parameters: 
- sgrnas: dataframe file of guides
- gene: fasta file of genes

### bigscam.sgrna.generate()
find_all and generate steps combined.
```
bigscam.sgrna.generate_BE(args*)
bigscam.sgrna.generate_DSB(args*)
bigscam.sgrna.generate_CRISPOR(args*)
```
Parameters: same as above


## Analysis

### bigscam.analysis.count_reads()
QC function to make sure all guides are in the library and quantify the skew present. Metrics are outputed. 
```
bigscam.analysis.count_reads(fastq, sgrnas)
```
Parameters: 
- fastq: data from NGS run
- sgrnas: dataframe of all sgRNAs

### bigscam.analysis.score()
```
bigscam.analysis.score()
```

## Plot

### bigscam.plot.coverage()
```
bigscam.plot.coverage()
```
Parameters: 
- a

### bigscam.plot.scattter()
```
bigscam.plot.scatter(counts)
```
Parameters: 
- a
![scatter1](https://github.com/liaulab/bigscam/assets/68132984/e46c2d96-fea4-4a1a-a217-3739168a9f79)
![scatter2](https://github.com/liaulab/bigscam/assets/68132984/ca2c33f0-073c-4a1d-be64-d43b669c7305)


### bigscam.plot.enrichment_map()
```
```
Parameters: 
- a
![enrichment_map](https://github.com/liaulab/bigscam/assets/68132984/32a92334-c749-476d-8dd4-d87f5e906456)

### bigscam.plot.volcano()
```
```
Parameters: 
- a

### bigscam.plot.pymol()
```
bigscam.plot.pymol(counts)
```
Parameters: 
- a


## Control Functions

### bigscam.control.print()
