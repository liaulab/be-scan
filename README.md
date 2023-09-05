# BIGSCAM (tentatively)
Base-editing Indels and Genomic Screen with CRISPR And More \
BioInformatics, Genomic Screening, and COmputational Analysis for Mutational data
Authors: SPS, XYH \

Description: We are building a package for all of the lab's computational needs. ^_^ This package is compatible with CRISPOR generated data.



# Documentation



## sgrna

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
Queries the CRISPOR server for a list of guides. 

### bigscam.sgrna.annotate()
```
bigscam.sgrna.annotate_BE(sgrnas, gene)
bigscam.sgrna.annotate_DSB(sgrnas, gene)
bigscam.sgrna.annotate_CRISPOR(sgrnas, gene)
```
### bigscam.sgrna.generate()

```
bigscam.sgrna.generate_BE(args*)
bigscam.sgrna.generate_DSB(args*)
bigscam.sgrna.generate_CRISPOR(args*)
```

```
bigscam.analysis.count_reads(fastq)
```

```
bigscam.analysis.score(counts)
```

```
bigscam.plot.scatter(counts)
```

```
bigscam.plot.pymol(counts)
```

```
```

```
```
