"""
Builders for synthetic gene .fasta files in the format be_scan expects.

The point of these is control. A real gene gives you whatever junction phases
it happens to have -- FOXA1, for instance, has two exons and therefore exactly
one internal junction. To exercise the intron/exon junction math properly we
need genes whose exon boundaries fall at every codon phase, with codons
deliberately split 1-2 and 2-1 across the boundary.

Input format, as parsed by GeneForCRISPR:

    >name range=chrN:START-END 5'pad=P 3'pad=P strand=+ repeatMasking=none
    <lowercase 5' pad><UPPERCASE EXON><lowercase 3' pad>

  * one record per exon
  * lowercase means intron/UTR, uppercase means coding
  * GeneForCRISPR.parse_exons reads the pad length from the first record and
    assumes it is the same everywhere
  * extract_metadata reads the chromosome from field 1 and the strand from
    field 4 of the first header line
"""

from be_scan.sgrna._genomic_ import DNA_AA_map

# Codons chosen so that single A->G or C->T edits produce clearly identifiable
# amino acid changes, and so the sequence contains no accidental stop codons.
SENSE_CODONS = [
    "ATG", "CTG", "GAA", "AAG", "ACC", "TGC", "GAT", "CAA", "GAG", "TAC",
    "AAA", "CCA", "ATC", "GCA", "TCA", "CAC", "AAC", "GTA", "ACA", "GGA",
    "CAT", "ATA", "TTA", "CAG", "AAT", "GCC", "TTG", "CCT", "ACT", "GTT",
    "CGA", "AGA", "GGC", "TGG", "CTA", "ATT", "CCG", "AGC", "TCC", "GTC",
]

STOP = "TAA"


def build_cds(n_codons=40):
    """
    Build a coding sequence that starts with ATG and ends with a stop codon.

    Returns (cds, protein) where protein includes the terminal '.' stop, so
    protein[i] is residue i+1 in be_scan's 1-indexed convention.
    """
    codons = [SENSE_CODONS[i % len(SENSE_CODONS)] for i in range(n_codons - 1)]
    codons[0] = "ATG"
    cds = "".join(codons) + STOP
    protein = "".join(DNA_AA_map[cds[i:i + 3]] for i in range(0, len(cds), 3))
    assert "." not in protein[:-1], "synthetic CDS has an internal stop codon"
    return cds, protein


def split_cds(cds, exon_lengths):
    """Cut a CDS into exons of the given lengths. Must sum to len(cds)."""
    assert sum(exon_lengths) == len(cds), (
        f"exon lengths sum to {sum(exon_lengths)} but CDS is {len(cds)} long"
    )
    exons, at = [], 0
    for length in exon_lengths:
        exons.append(cds[at:at + length])
        at += length
    return exons


def _intron(length, donor=True, acceptor=True, filler="c"):
    """
    A lowercase intron fragment.

    Real introns begin 'gt' and end 'ag'; be_scan's splice-site detection keys
    off exactly those two dinucleotides, so they must be present for
    is_splice_donor_guide / is_splice_acceptor_guide to be exercised.
    """
    body = filler * length
    seq = list(body)
    if donor:
        seq[0], seq[1] = "g", "t"
    if acceptor:
        seq[-2], seq[-1] = "a", "g"
    return "".join(seq)


def build_fasta(
    exons,
    pad=30,
    strand="+",
    chrom="chr1",
    start=100000,
    name="synthetic",
    utr_5="acgtacgtacgtacgtacgtacgtacgtac",
    utr_3="tgcatgcatgcatgcatgcatgcatgcatg",
):
    """
    Assemble a .fasta string, one record per exon.

    The 5' pad of the first exon and the 3' pad of the last exon are UTR
    sequence. Every internal pad is intron: the 3' pad of exon e opens with
    'gt' (splice donor) and the 5' pad of exon e+1 closes with 'ag' (splice
    acceptor).
    """
    assert len(utr_5) == pad and len(utr_3) == pad, "UTR pads must be `pad` long"

    records, position = [], start
    for e, exon in enumerate(exons):
        five = utr_5 if e == 0 else _intron(pad, donor=False, acceptor=True)
        three = utr_3 if e == len(exons) - 1 else _intron(pad, donor=True, acceptor=False)
        body = five + exon.upper() + three

        end = position + len(body) - 1
        header = (
            f">{name}_{e} range={chrom}:{position}-{end} "
            f"5'pad={pad} 3'pad={pad} strand={strand} repeatMasking=none"
        )
        # wrap at 50 columns, as UCSC does
        wrapped = "\n".join(body[i:i + 50] for i in range(0, len(body), 50))
        records.append(f"{header}\n{wrapped}")
        position = end + 1

    return "\n".join(records) + "\n"


def write_gene(tmpdir, exon_lengths, n_codons=40, pad=30, strand="+", name="synthetic"):
    """
    Build a complete synthetic gene and write it to tmpdir.

    Returns a dict with everything a test needs to check annotations against:
      path        -- the written .fasta
      cds         -- the coding sequence
      protein     -- translation, 1-indexed as protein[pos-1]
      exons       -- list of exon sequences (uppercase)
      pad         -- intron/UTR pad length
      junctions   -- CDS offset of each internal exon boundary
      phases      -- offset % 3 for each internal boundary
    """
    cds, protein = build_cds(n_codons)
    exons = split_cds(cds, exon_lengths)
    fasta = build_fasta(exons, pad=pad, strand=strand, name=name)

    path = tmpdir / f"{name}.fasta"
    path.write_text(fasta)

    offsets, at = [], 0
    for length in exon_lengths[:-1]:
        at += length
        offsets.append(at)

    return {
        "path": path,
        "fasta": fasta,
        "cds": cds,
        "protein": protein,
        "exons": exons,
        "exon_lengths": list(exon_lengths),
        "pad": pad,
        "junctions": offsets,
        "phases": [o % 3 for o in offsets],
    }


# Exon length sets that place internal boundaries at every codon phase.
# Cumulative offsets and their phases are given in the comments.
ALL_PHASE_LENGTHS = [30, 31, 34, 25]  # 30->0, 61->1, 95->2
PHASE0_LENGTHS = [30, 30, 60]         # 30->0, 60->0
PHASE1_LENGTHS = [31, 30, 59]         # 31->1, 61->1
PHASE2_LENGTHS = [32, 30, 58]         # 32->2, 62->2
SINGLE_EXON_LENGTHS = [120]           # no internal junction, control case
