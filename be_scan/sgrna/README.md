# `be_scan.sgrna` — library design and annotation

Turns a gene sequence into an annotated sgRNA library: which guides exist,
and which residues each one could change.

## Public functions

| function | purpose |
| --- | --- |
| `generate_library` | find every guide with a matching PAM, optionally requiring an editable base |
| `annotate` | predict the amino acid changes each guide could produce |
| `reference_check` | count occurrences of each spacer in a reference genome |
| `coverage_plots` | summarise how much of the protein the library reaches |
| `design_library` | run all four in order |

## Input format

A UCSC-style `.fasta`, one record per exon, exons uppercase and flanking
intron lowercase:

```
>NAME_0 range=chr14:37594871-37595002 5'pad=30 3'pad=30 strand=- repeatMasking=none
gatggttgtattgggcagggtggctccaggATGTTAGGAACTGTGAAGAT
...AGgtgagaggcaaggcgaggcgggcggccggt
```

The case distinction is load-bearing: it is how the package knows where the
coding sequence starts and stops. Padding must be the same length on both
sides of every exon.

## Quick start

```python
from be_scan.sgrna import generate_library, annotate

guides, gene = generate_library(
    gene_filepath="FOXA1.fasta", gene_name="FOXA1",
    cas_type="SpG",                 # NGN
    edit_from="A", edit_to="G",     # ABE
    window=[4, 8],
)
annotated = annotate(guides_file=guides, edit_from="A", edit_to="G",
                     exons=gene.exons, window=[4, 8])
```

Annotating repeatedly adds a prefixed block per editor, so `A`→`G`, `C`→`T`
and the dual `AC`→`GT` can share one table as `AtoG_*`, `CtoT_*`, `ACtoGT_*`.

## Output columns

Per editor prefix:

| column | meaning |
| --- | --- |
| `*_win_overlap` | `Exon`, `Intron` or `Exon/Intron` |
| `*_mutations` | predicted changes; `;` separates alternative outcomes, `/` separates residues changed together |
| `*_muttypes` | the class of each call above |
| `*_muttype` | one label, by priority: splice site, nonsense, missense, silent, then non-coding |
| `*_pos` | residue position, or `-1` where not applicable |
| `<base>_count` | editable bases in the window |

## `window`

The editing window in 1-indexed guide coordinates, `[4, 8]` by default.

Any window with **`start >= 3` and `end <= 18`** is supported. Outside that
range a codon overlapping guide positions 1–2 or 19–20 can extend past the end
of the guide and cannot be resolved from the guide sequence alone; those bases
are skipped rather than mis-numbered.

Before 3.0.0 only `[4, 8]` was correct — see the changelog.

## `guide_len`

Defaults to 20. The searched sequence is `guide_len + len(PAM)`, so shorter
spacers and non-3 nt PAMs both work.

## Intron/exon junctions

A codon split across a junction is completed from the neighbouring exon before
translation. This is the part of the module most worth trusting carefully, so
it is cross-checked in the test suite against a second implementation that
shares none of its arithmetic — one works in guide-relative frames, the other
in absolute CDS coordinates. They agree on 33,801 mutation calls across three
editors on real input, and on synthetic genes built with junctions at all
three codon phases.

If you supply `protein_filepath`, every predicted wild-type residue is checked
against that sequence and a mismatch warns. Worth doing once per gene: a frame
error shows up immediately.

## Custom column names

`annotate` takes `seq_col`, `frame_col`, `strand_col`, `gene_pos_col`,
`window_start_col`, `window_end_col` and `gene_col`. Before 3.0.0 several of
these were accepted and then ignored.
