"""
Cross-check be_scan's annotation math against an independent implementation.

Every guide produced by generate_library is annotated twice: once by
be_scan.sgrna.annotate (frame arithmetic, codon-window slicing, junction
patch-up) and once by tests.sgrna_test.reference_annotator (absolute CDS
coordinates, whole-CDS retranslation). The two must agree.

These tests do not modify baseline behavior. Where they disagree with the
baseline, they report the disagreement.
"""

import warnings

import pandas as pd
import pytest

from be_scan.sgrna._gene_ import GeneForCRISPR
from be_scan.sgrna._genomic_ import complements, rev_complement
from be_scan.sgrna.annotate import annotate
from be_scan.sgrna.generate_library import generate_library

from tests.sgrna_test import reference_annotator as ref
from tests.sgrna_test import synthetic_gene as syn

N_MER = 23
GUIDE_LEN = 20

EDITORS = [
    pytest.param("A", "G", id="ABE_AtoG"),
    pytest.param("C", "T", id="CBE_CtoT"),
    pytest.param("AC", "GT", id="dual_ACtoGT"),
]

GENE_LAYOUTS = [
    pytest.param(syn.ALL_PHASE_LENGTHS, id="junctions_phase_0_1_2"),
    pytest.param(syn.PHASE0_LENGTHS, id="junctions_phase_0"),
    pytest.param(syn.PHASE1_LENGTHS, id="junctions_phase_1"),
    pytest.param(syn.PHASE2_LENGTHS, id="junctions_phase_2"),
    pytest.param(syn.SINGLE_EXON_LENGTHS, id="single_exon_no_junction"),
]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def guide_positions(gene):
    """
    Recreate the (record, index, strand) -> spacer enumeration that
    GeneForCRISPR.find_all_guides performs, so an annotated row can be traced
    back to the gene coordinates it came from.

    Sequences that occur more than once are dropped: generate_library removes
    duplicate spacers by default, so an ambiguous sequence cannot be matched
    back to a single origin.
    """
    seen = {}
    duplicated = set()
    for e, record in enumerate(gene.exons_extra):
        for i in range(len(record) - N_MER - 1):
            block = record[i:i + N_MER]
            fwd = block[:GUIDE_LEN]
            rev = rev_complement(complements, block[3:])
            for seq, strand in ((fwd, "sense"), (rev, "antisense")):
                if seq in seen and seen[seq] != (e, i, strand):
                    duplicated.add(seq)
                seen.setdefault(seq, (e, i, strand))
    for seq in duplicated:
        seen.pop(seq, None)
    return seen


def build(tmp_path, exon_lengths, edit_from, edit_to, window=(4, 8),
          cas_type="SpRY", exclude_sequences=()):
    """Build a synthetic gene, generate its library, and annotate it."""
    info = syn.write_gene(tmp_path, exon_lengths)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        guides, gene = generate_library(
            gene_filepath=str(info["path"]),
            gene_name="SYN",
            cas_type=cas_type,
            edit_from=edit_from,
            edit_to=edit_to,
            window=list(window),
            return_df=True,
            save_df=False,
            exclude_sequences=list(exclude_sequences),
        )
        annotated = annotate(
            guides_file=guides,
            edit_from=edit_from,
            edit_to=edit_to,
            exons=gene.exons,
            window=list(window),
            return_df=True,
            save_df=False,
        )

    info["gene"] = gene
    info["guides"] = guides
    info["annotated"] = annotated
    info["lookup"] = guide_positions(gene)
    info["coord_maps"] = ref.build_coordinate_map(info["exons"], info["pad"])
    return info


def compare_all(info, edit_from, edit_to, window=(4, 8)):
    """
    Compare be_scan's mutation calls against the reference for every guide.

    Returns (n_compared, disagreements) where each disagreement records the
    guide and both sets of predicted mutations.
    """
    pre = f"{edit_from}to{edit_to}"
    disagreements = []
    compared = 0

    for _, row in info["annotated"].iterrows():
        seq = row["sgRNA_seq"]
        origin = info["lookup"].get(seq)
        if origin is None:
            continue
        e, i, strand = origin
        assert strand == row["sgRNA_strand"], (
            f"strand mismatch for {seq}: enumeration says {strand}, "
            f"be_scan says {row['sgRNA_strand']}"
        )

        expected = ref.expected_mutations(
            record=info["gene"].exons_extra[e],
            index=i,
            window=window,
            edit_from=edit_from,
            edit_to=edit_to,
            strand=strand,
            coord_map=info["coord_maps"][e],
            cds=info["cds"],
            protein=info["protein"],
        )
        actual = ref.parse_mutations(row[f"{pre}_mutations"])
        compared += 1

        if expected != actual:
            disagreements.append({
                "sgRNA_ID": row["sgRNA_ID"],
                "seq": seq,
                "exon": e,
                "index": i,
                "strand": strand,
                "frame": row["starting_frame"],
                "gene_pos": row["gene_pos"],
                "win_start": row["windowstart_pos"],
                "win_end": row["windowend_pos"],
                "win_overlap": row[f"{pre}_win_overlap"],
                "expected": sorted("/".join(sorted(c)) for c in expected),
                "actual": sorted("/".join(sorted(c)) for c in actual),
            })

    return compared, disagreements


def format_disagreements(disagreements, limit=10):
    lines = []
    for d in disagreements[:limit]:
        lines.append(
            f"  {d['sgRNA_ID']} {d['seq']} exon={d['exon']} idx={d['index']} "
            f"{d['strand']} frame={d['frame']} overlap={d['win_overlap']}\n"
            f"      reference: {d['expected']}\n"
            f"      be_scan  : {d['actual']}"
        )
    if len(disagreements) > limit:
        lines.append(f"  ... and {len(disagreements) - limit} more")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# the synthetic gene fixtures themselves must be sane
# --------------------------------------------------------------------------- #


def test_synthetic_cds_is_well_formed(tmp_path):
    info = syn.write_gene(tmp_path, syn.ALL_PHASE_LENGTHS)
    assert len(info["cds"]) % 3 == 0
    assert info["cds"].startswith("ATG")
    assert info["protein"].endswith(".")
    assert "." not in info["protein"][:-1], "internal stop in synthetic CDS"
    assert info["phases"] == [0, 1, 2], (
        f"layout should place junctions at all three phases, got {info['phases']}"
    )


def test_parsed_exons_match_the_source(tmp_path):
    """GeneForCRISPR must recover exactly the exons we wrote."""
    info = syn.write_gene(tmp_path, syn.ALL_PHASE_LENGTHS)
    gene = GeneForCRISPR(filepath=str(info["path"]))
    gene.parse_exons()
    assert gene.exons == info["exons"]
    assert "".join(gene.exons) == info["cds"]
    assert gene.intron_len == info["pad"]


# --------------------------------------------------------------------------- #
# the main cross-check
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("exon_lengths", GENE_LAYOUTS)
@pytest.mark.parametrize("edit_from,edit_to", EDITORS)
def test_mutations_match_reference(tmp_path, exon_lengths, edit_from, edit_to):
    info = build(tmp_path, exon_lengths, edit_from, edit_to)
    compared, disagreements = compare_all(info, edit_from, edit_to)

    assert compared > 50, f"only {compared} guides compared; fixture too small"
    assert not disagreements, (
        f"{len(disagreements)}/{compared} guides disagree with the reference "
        f"implementation:\n{format_disagreements(disagreements)}"
    )


@pytest.mark.parametrize(
    "exon_lengths",
    [
        pytest.param([30, 31, 34, 25], id="layout_A"),
        pytest.param([27, 35, 33, 25], id="layout_B"),
        pytest.param([21, 22, 23, 54], id="layout_C"),
        pytest.param([45, 20, 20, 35], id="layout_D"),
        pytest.param([60, 29, 31], id="layout_E"),
        pytest.param([24, 24, 24, 24, 24], id="layout_F_five_exons"),
    ],
)
def test_default_window_holds_across_gene_layouts(tmp_path, exon_lengths):
    """
    The default window [4,8] must be correct for any exon layout.

    Different layouts place the junctions at different codon phases and split
    codons 1-2 or 2-1 across the boundary, which is what fill_in_dna_exons has
    to reassemble.
    """
    info = build(tmp_path, exon_lengths, "A", "G")
    compared, disagreements = compare_all(info, "A", "G")

    assert compared > 50, f"only {compared} guides compared"
    assert not disagreements, (
        f"layout {exon_lengths}: {len(disagreements)}/{compared} guides "
        f"disagree:\n{format_disagreements(disagreements)}"
    )


# --------------------------------------------------------------------------- #
# characterization: the window parameter only works at its default
# --------------------------------------------------------------------------- #

# `window` is documented as tunable in generate_library, annotate and
# design_library ("Editing window, 4th to 8th bases inclusive by default").
# It is not. Sweeping every window [start, start+width-1] for start 1..10 and
# width 1..7 against the reference implementation, exactly one of the 70
# combinations is correct, and it is the default:
#
#   start |    w=1    w=2    w=3    w=4    w=5    w=6    w=7
#       1 |  CRASH  CRASH  CRASH  CRASH  CRASH  CRASH  CRASH
#       2 |  CRASH  CRASH  CRASH  CRASH  CRASH  CRASH  CRASH
#       3 |  CRASH  CRASH  CRASH  CRASH  CRASH  CRASH  CRASH
#       4 |  CRASH  CRASH  CRASH  CRASH     OK  WRONG  WRONG
#       5 |  CRASH  CRASH  CRASH  CRASH  WRONG  WRONG  CRASH
#       6 |  CRASH  CRASH  CRASH  CRASH  WRONG  CRASH  CRASH
#    7-10 |  CRASH  CRASH  CRASH  CRASH  CRASH  CRASH  CRASH
#
# Root causes, all in code tuned to window == [4, 8]:
#
#   _gene_.py:88,121   `seq[window[1]]` indexes one past the window's last
#                      base. On a 20-mer reverse guide this raises IndexError
#                      whenever window[1] == 20.
#   _guideRNA_.py:117  calc_target's codon window is anchored to the guide
#                      start, not to the editing window, so a window far from
#                      the 5' end falls outside the translated codons and
#                      DNA_to_AA asserts on a length not divisible by 3.
#   _guideRNA_.py:161  the residue-index arithmetic in annotate_muts hardcodes
#                      offsets (+3, <3, +1) calibrated for [4,8]; other
#                      windows shift the reported residue by one.
#   _guideRNA_.py:255  fill_in_dna_exons reaches exons[exon_i +/- 1] without a
#                      bounds check, so a boundary codon in the first or last
#                      exon indexes off the end of the list.
#
# These are recorded, not fixed: the baseline annotation logic is deliberately
# left untouched. If any of these xfails starts passing, the underlying issue
# was addressed and this block should be updated.

WINDOW_CRASHES = [(1, 5), (2, 6), (3, 7), (3, 9), (5, 7), (4, 4), (8, 12), (1, 20)]
WINDOW_WRONG = [(5, 9), (6, 10), (4, 9), (4, 10), (4, 11)]


@pytest.mark.parametrize("window", WINDOW_CRASHES, ids=str)
@pytest.mark.xfail(reason="baseline: non-default editing windows raise", strict=True)
def test_nondefault_window_does_not_crash(tmp_path, window):
    """Characterizes that non-default windows raise rather than annotate."""
    build(tmp_path, syn.ALL_PHASE_LENGTHS, "A", "G", window=window)


@pytest.mark.parametrize("window", WINDOW_WRONG, ids=str)
@pytest.mark.xfail(
    reason="baseline: non-default windows shift residue numbering", strict=True
)
def test_nondefault_window_matches_reference(tmp_path, window):
    """Characterizes that non-default windows misnumber residues."""
    info = build(tmp_path, syn.ALL_PHASE_LENGTHS, "A", "G", window=window)
    _, disagreements = compare_all(info, "A", "G", window=window)
    assert not disagreements, format_disagreements(disagreements)


def test_default_window_is_the_only_supported_window():
    """
    A signpost, so the [4,8]-only limitation is discoverable from the tests.

    Kept deliberately simple: the detailed evidence is in the sweep recorded
    in the comment block above.
    """
    from be_scan.sgrna.annotate import annotate as _annotate
    import inspect

    default = inspect.signature(_annotate).parameters["window"].default
    assert default == [4, 8], (
        "annotate's default window changed. The annotation math is only "
        "correct at [4,8]; see the window sweep above."
    )


@pytest.mark.parametrize("edit_from,edit_to", EDITORS)
def test_junction_spanning_guides_match_reference(tmp_path, edit_from, edit_to):
    """
    Restrict the comparison to guides whose window straddles a junction.

    These are the guides that exercise fill_in_dna_exons, where a codon is
    completed using bases from the neighbouring exon.
    """
    info = build(tmp_path, syn.ALL_PHASE_LENGTHS, edit_from, edit_to)
    pre = f"{edit_from}to{edit_to}"

    junction = info["annotated"][
        info["annotated"][f"{pre}_win_overlap"] == "Exon/Intron"
    ]
    assert len(junction) > 0, "fixture produced no junction-spanning guides"

    _, disagreements = compare_all(info, edit_from, edit_to)
    junction_seqs = set(junction["sgRNA_seq"])
    junction_bad = [d for d in disagreements if d["seq"] in junction_seqs]

    assert not junction_bad, (
        f"{len(junction_bad)}/{len(junction)} junction-spanning guides "
        f"disagree with the reference:\n{format_disagreements(junction_bad)}"
    )
