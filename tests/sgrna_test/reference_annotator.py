"""
An independent reimplementation of base-editor mutation annotation.

be_scan derives mutations through frame arithmetic: it tracks a starting_frame
per guide, slices a codon window out of the 23-mer, and patches codons that
straddle an intron/exon boundary back together from the neighbouring exon
(_guideRNA_.fill_in_dna_exons). That arithmetic is exactly what these tests
are meant to check, so a reference implementation that reused any of it would
be worthless.

This module takes the opposite route. It never computes a frame. It works in
gene coordinates, maps each edited base to an absolute CDS index, rebuilds the
whole coding sequence, translates it, and diffs against the wild-type protein.
A codon split across an intron/exon junction needs no special handling here,
because in CDS space it was never split.

Where the two approaches agree, the frame arithmetic is right. Where they
disagree, the disagreement is reported -- the baseline code is not modified.
"""

from itertools import product

from be_scan.sgrna._genomic_ import DNA_AA_map, complements

GUIDE_LEN = 20
N_MER = 23  # spacer + PAM, as used by GeneForCRISPR.find_all_guides


def translate(seq):
    """Translate a coding sequence. Assumes length is a multiple of 3."""
    return "".join(DNA_AA_map[seq[i:i + 3].upper()] for i in range(0, len(seq), 3))


def build_coordinate_map(exons, pad):
    """
    Map (record index, position within record) -> absolute CDS index.

    Each fasta record is <pad lowercase><EXON uppercase><pad lowercase>, and
    guides never span two records. A position is coding only if it falls in
    the uppercase middle.

    Returns a list of dicts, one per record, mapping local position -> CDS index.
    """
    maps, cds_offset = [], 0
    for exon in exons:
        maps.append({pad + k: cds_offset + k for k in range(len(exon))})
        cds_offset += len(exon)
    return maps


def window_positions(index, window, strand):
    """
    Absolute positions within a record covered by the editing window.

    A sense guide's spacer occupies record[index : index+20], so window
    position w (1-indexed) is record[index + w - 1].

    An antisense guide is built by GeneForCRISPR as
    revcomp(record[index+3 : index+23]), so its spacer position k (0-indexed)
    corresponds to record position index + 22 - k.
    """
    if strand == "sense":
        return [index + w - 1 for w in range(window[0], window[1] + 1)]
    return [index + (N_MER - 1) - (w - 1) for w in range(window[0], window[1] + 1)]


def editable_positions(record, positions, edit_from, strand):
    """
    Which of these record positions can the editor actually change?

    be_scan matches edit bases case-sensitively against uppercase letters, so
    lowercase (intronic) bases are never edited. An antisense A->G editor
    changes a sense-strand T to a C, so on the sense strand we look for the
    complement of edit_from.
    """
    if strand == "sense":
        targets = set(edit_from.upper())
    else:
        targets = {complements[b] for b in edit_from.upper()}
    return [p for p in positions if record[p] in targets]


def _replacement_map(edit_from, edit_to, strand):
    """Sense-strand base -> the set of bases it may become."""
    mapping = {}
    for src, dst in zip(edit_from.upper(), edit_to.upper()):
        if strand == "sense":
            a, b = src, dst
        else:
            a, b = complements[src], complements[dst]
        mapping.setdefault(a, set()).update({a, b})
    return mapping


def expected_mutations(record, index, window, edit_from, edit_to, strand,
                       coord_map, cds, protein):
    """
    Every mutation combination this guide could produce.

    Returns a set of frozensets. Each frozenset is one edit combination,
    holding 'X###Y' strings for the residues that combination changes. The
    empty combination (edit nothing) is excluded, matching be_scan, which
    skips combinations where the DNA is unchanged.
    """
    positions = window_positions(index, window, strand)
    positions = [p for p in positions if 0 <= p < len(record)]
    targets = editable_positions(record, positions, edit_from, strand)

    # Only bases inside an exon can change the protein. Intronic edits are
    # real but produce no coding consequence, which is what be_scan reports.
    coding = [(p, coord_map[p]) for p in targets if p in coord_map]
    if not coding:
        return set()

    replacements = _replacement_map(edit_from, edit_to, strand)
    choices = [sorted(replacements[record[p]]) for p, _ in coding]

    results = set()
    for combo in product(*choices):
        new_cds = list(cds)
        changed_codons = set()
        for (pos, cds_index), base in zip(coding, combo):
            if base != cds[cds_index]:
                changed_codons.add(cds_index // 3)
            new_cds[cds_index] = base
        if not changed_codons:
            continue  # be_scan skips combinations that leave the DNA identical

        new_protein = translate("".join(new_cds))
        # be_scan's format_mutation emits a call for every codon whose DNA
        # changed, whether or not the residue changed. A silent edit is
        # reported as e.g. 'L2L', so codon identity -- not residue identity --
        # is the right comparison.
        muts = frozenset(
            f"{protein[c]}{c + 1}{new_protein[c]}" for c in sorted(changed_codons)
        )
        results.add(muts)
    return results


def parse_mutations(cell):
    """
    Parse be_scan's mutation string into a set of frozensets.

    Format is 'A1B/C2D;E3F': ';' separates alternative edit combinations, '/'
    separates residues within one combination. be_scan builds this with
    ';'.join(set(...)), so ordering is not stable and comparison must be
    order-insensitive.
    """
    if not isinstance(cell, str) or not cell.strip():
        return set()
    out = set()
    for combo in cell.split(";"):
        residues = frozenset(m for m in combo.split("/") if m)
        if residues:
            out.add(residues)
    return out


def record_strings(gene):
    """The per-exon record strings held by a parsed GeneForCRISPR."""
    return gene.exons_extra
