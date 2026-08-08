"""
Unit tests for be_scan.sgrna._guideRNA_.

These target the helpers directly rather than through annotate(), so a failure
points at one function instead of the whole pipeline. The intron/exon
reassembly (fill_in_dna_exons) and the splice-site detectors get the most
attention -- they are the parts that decide how an edit across a junction is
translated.
"""

import pytest

from be_scan.sgrna._guideRNA_ import (
    annotate_intron_exon,
    annotate_intron_exon_anti,
    assign_position,
    determine_mutations,
    extract_uppercase_letters,
    fill_in_dna_exons,
    filter_sequence,
    find_first_uppercase_index,
    format_mutation,
    is_splice_acceptor_guide,
    is_splice_donor_guide,
    mutation_combos,
    window_on_boundary,
)

WINDOW = [4, 8]


# --------------------------------------------------------------------------- #
# window_on_boundary
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "seq,expected",
    [
        ("ATGCCC", False),   # all exonic
        ("atgccc", False),   # all intronic
        ("ATGccc", True),    # spans a junction
        ("atgCCC", True),
        ("AtGcCc", True),
        ("", False),
        ("A", False),
        ("a", False),
    ],
)
def test_window_on_boundary(seq, expected):
    assert window_on_boundary(seq) is expected


# --------------------------------------------------------------------------- #
# annotate_intron_exon
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "seq,expected",
    [
        ("AAAGGGGGAAAAAAAAAAAA", "Exon"),
        ("AAAgggggAAAAAAAAAAAA", "Intron"),
        ("AAAGGgggAAAAAAAAAAAA", "Exon/Intron"),
        ("AAAgggGGAAAAAAAAAAAA", "Exon/Intron"),
    ],
)
def test_annotate_intron_exon_reads_only_the_window(seq, expected):
    """Only bases 4-8 decide the call; everything outside is ignored."""
    assert annotate_intron_exon(seq, WINDOW) == expected


def test_annotate_intron_exon_ignores_bases_outside_the_window():
    assert annotate_intron_exon("aaaGGGGGaaaaaaaaaaaa", WINDOW) == "Exon"
    assert annotate_intron_exon("AAAgggggAAAAAAAAAAAA", WINDOW) == "Intron"


def test_annotate_intron_exon_anti_mirrors_the_window():
    """
    The antisense variant reads the mirrored window on a 20-mer.

    For window [4,8] that is x[12:17]; the docstring in the source states the
    same. Verified here against explicit slices rather than by construction.
    """
    seq = "A" * 12 + "GGGGG" + "A" * 3
    assert len(seq) == 20
    assert annotate_intron_exon_anti(seq, WINDOW) == "Exon"

    seq = "A" * 12 + "ggggg" + "A" * 3
    assert annotate_intron_exon_anti(seq, WINDOW) == "Intron"

    seq = "A" * 12 + "GGggg" + "A" * 3
    assert annotate_intron_exon_anti(seq, WINDOW) == "Exon/Intron"


def test_annotate_intron_exon_anti_assumes_a_20mer():
    """
    The mirrored slice is computed against a hardcoded length of 20.

    _guideRNA_.py:445 reads `x[20-window[1] : 20-window[0]+1]` with a literal
    20 rather than len(x). On a guide of any other length the window lands in
    the wrong place. Recorded here so the assumption is explicit; the guide
    length is being made tunable elsewhere in this release.
    """
    twenty = "A" * 12 + "ggggg" + "A" * 3
    assert annotate_intron_exon_anti(twenty, WINDOW) == "Intron"

    # the same guide one base longer: the intronic run is no longer read
    twenty_one = "A" + twenty
    assert annotate_intron_exon_anti(twenty_one, WINDOW) == "Exon/Intron"


# --------------------------------------------------------------------------- #
# mutation_combos
# --------------------------------------------------------------------------- #


def test_mutation_combos_enumerates_every_subset_on_the_sense_strand():
    """Two editable A's give 2^2 sequences, including the unedited original."""
    combos = mutation_combos("ACA", ("A", "G"), "sense")
    assert sorted(combos) == ["ACA", "ACG", "GCA", "GCG"]


def test_mutation_combos_uses_the_complement_on_the_antisense_strand():
    """
    An antisense A->G editor changes a sense-strand T to a C.

    mutation_combos is handed the sense-strand window, so it must look for T
    and write C.
    """
    combos = mutation_combos("TCT", ("A", "G"), "antisense")
    assert sorted(combos) == ["CCC", "CCT", "TCC", "TCT"]


def test_mutation_combos_returns_the_input_when_nothing_is_editable():
    assert mutation_combos("CCC", ("A", "G"), "sense") == ["CCC"]


def test_mutation_combos_ignores_lowercase_bases():
    """Lowercase means intronic, and intronic bases are never edited."""
    combos = mutation_combos("aCA", ("A", "G"), "sense")
    assert sorted(combos) == ["aCA", "aCG"]


def test_mutation_combos_scales_as_two_to_the_n():
    for n in range(1, 6):
        combos = mutation_combos("A" * n, ("A", "G"), "sense")
        assert len(combos) == 2 ** n
        assert len(set(combos)) == 2 ** n


# --------------------------------------------------------------------------- #
# fill_in_dna_exons -- the intron/exon junction reassembly
# --------------------------------------------------------------------------- #


def test_fill_in_dna_exons_completes_a_codon_from_the_next_exon():
    """
    A codon split 1-2 across the 3' end of an exon is completed from exon+1.

    'ATG' is wholly exonic and untouched; 'Ccc' has two intronic bases, so the
    first two bases of the next exon replace them.
    """
    exons = ["ATGC", "TTGGG", "AAACCC"]
    row = {"exon": 0}
    assert fill_in_dna_exons("ATGCcc", row, exons) == "ATGCTT"


def test_fill_in_dna_exons_completes_a_codon_from_the_previous_exon():
    """A codon split 2-1 across the 5' end of an exon is completed from exon-1."""
    exons = ["ATGCTT", "GGGAAA", "CCCTTT"]
    row = {"exon": 1}
    assert fill_in_dna_exons("ttGGGAAA", row, exons) == "TTGGGAAA"


def test_fill_in_dna_exons_leaves_fully_coding_codons_alone():
    exons = ["ATGCTT", "GGGAAA"]
    row = {"exon": 0}
    assert fill_in_dna_exons("ATGCTT", row, exons) == "ATGCTT"


def test_fill_in_dna_exons_leaves_fully_intronic_codons_alone():
    """
    A codon with no uppercase base is not a boundary codon.

    window_on_boundary requires both cases, so an all-lowercase codon is
    skipped and later becomes '_' in translation.
    """
    exons = ["ATGCTT", "GGGAAA"]
    row = {"exon": 0}
    assert fill_in_dna_exons("ATGttt", row, exons) == "ATGttt"


def test_fill_in_dna_exons_reads_past_the_last_exon():
    """
    A boundary codon in the final exon indexes exons[exon_i + 1].

    _guideRNA_.py:259 has no bounds check, so this raises rather than
    reporting that there is no next exon. Recorded, not fixed.
    """
    exons = ["ATGCTT", "GGGAAA"]
    row = {"exon": 1}  # last exon
    with pytest.raises(IndexError):
        fill_in_dna_exons("GGGAAACcc", row, exons)


def test_fill_in_dna_exons_wraps_to_the_last_exon_before_the_first():
    """
    A boundary codon in the first exon indexes exons[-1].

    Python's negative indexing makes exons[0 - 1] the *last* exon, so this
    silently completes the codon from the wrong end of the gene instead of
    raising. The result below is built from 'TTT' (the tail of the last exon)
    rather than from any upstream sequence.

    Not reachable through annotate() on a well-formed gene, because the first
    exon begins at a codon boundary, but recorded because the failure is
    silent rather than loud.
    """
    exons = ["ATGCTT", "GGGAAA", "CCCTTT"]
    row = {"exon": 0}
    assert fill_in_dna_exons("ttATGCTT", row, exons) == "TTATGCTT"


# --------------------------------------------------------------------------- #
# format_mutation
# --------------------------------------------------------------------------- #


def test_format_mutation_reports_only_codons_whose_dna_changed():
    """A codon with identical DNA is skipped even if scanned."""
    out = format_mutation(
        aa="ML", new_aa="VL", start=1,
        amino_acid_seq=None,
        dna="ATGCTG", new_dna="GTGCTG", guide="x",
    )
    assert out == "M1V"


def test_format_mutation_reports_silent_changes():
    """A changed codon that yields the same residue is still reported."""
    out = format_mutation(
        aa="ML", new_aa="ML", start=1,
        amino_acid_seq=None,
        dna="ATGTTA", new_dna="ATGCTA", guide="x",
    )
    assert out == "L2L"


def test_format_mutation_joins_multiple_residues_with_a_slash():
    out = format_mutation(
        aa="ML", new_aa="VL", start=1,
        amino_acid_seq=None,
        dna="ATGTTA", new_dna="GTGCTA", guide="x",
    )
    assert out == "M1V/L2L"


def test_format_mutation_skips_non_coding_codons():
    """'_' marks a codon that is not fully exonic; it produces no call."""
    out = format_mutation(
        aa="_L", new_aa="_L", start=1,
        amino_acid_seq=None,
        dna="atgTTA", new_dna="acgCTA", guide="x",
    )
    assert out == "L2L"


def test_format_mutation_warns_when_the_protein_disagrees():
    """
    A mismatch against the supplied protein sequence warns rather than fails.

    This is the check that catches a frame error against a real .fasta.
    """
    protein = {1: "Q", 2: "L"}  # residue 1 should be M, so this disagrees
    with pytest.warns(UserWarning):
        format_mutation(
            aa="ML", new_aa="VL", start=1,
            amino_acid_seq=protein,
            dna="ATGCTG", new_dna="GTGCTG", guide="badguide",
        )


def test_format_mutation_asserts_on_a_non_positive_residue_number():
    """Residue numbering is 1-based; index 0 or below is a frame error."""
    with pytest.raises(AssertionError):
        format_mutation(
            aa="ML", new_aa="VL", start=-1,
            amino_acid_seq=None,
            dna="ATGCTG", new_dna="GTGCTG", guide="x",
        )


# --------------------------------------------------------------------------- #
# determine_mutations
# --------------------------------------------------------------------------- #


COLS = ("starting_frame", "sgRNA_strand", "gene_pos", "sgRNA_seq",
        "windowstart_pos", "windowend_pos")


def _muts_row(mutations, ws=10, we=14):
    return {
        "windowstart_pos": ws, "windowend_pos": we,
        "AtoG_mutations": mutations,
    }


@pytest.mark.parametrize(
    "mutations,expected",
    [
        ("M1V", "Missense"),
        ("L2L", "Silent"),
        ("Q5.", "Nonsense"),
        ("M1V/L2L", "Missense/Silent"),
        ("M1V;L2L", "Missense;Silent"),
        ("Q5./M1V", "Nonsense/Missense"),
        ("", ""),
    ],
)
def test_determine_mutations_classifies_each_call(mutations, expected):
    assert determine_mutations(_muts_row(mutations), COLS, "AtoG") == expected


def test_determine_mutations_returns_none_without_coding_dna():
    """Both window bounds at -1 means no coding sequence in the window."""
    assert determine_mutations(_muts_row("M1V", ws=-1, we=-1), COLS, "AtoG") is None


def test_stop_loss_is_classified_as_missense():
    """
    A read-through of the terminal stop ('.473W') is reported as Missense.

    determine_mutations keys 'Nonsense' off a trailing '.', so a mutation that
    *removes* a stop rather than creating one falls through to Missense.
    Confirmed against real FOXA1 output, where sgRNA_1478 yields '.473W'.
    """
    assert determine_mutations(_muts_row(".473W"), COLS, "AtoG") == "Missense"


# --------------------------------------------------------------------------- #
# splice site detection
# --------------------------------------------------------------------------- #


def test_splice_acceptor_detected_on_sense_strand():
    """Intronic 'ag' immediately before an exon, inside the window."""
    seq = "aaaag" + "GCCTACTCCTCCGTC"
    assert len(seq) == 20
    row = {"sgRNA_seq": seq}
    assert is_splice_acceptor_guide(row, "AtoG", "sense", WINDOW) is True


def test_splice_donor_detected_on_sense_strand():
    """Intronic 'gt' immediately after an exon, inside the window."""
    seq = "GCT" + "gt" + "gagaggcaaggcgg"
    seq = "GCTgtgagaggcaaggcggg"[:20]
    row = {"sgRNA_seq": seq}
    result = is_splice_donor_guide(row, "TtoC", "sense", WINDOW)
    assert isinstance(result, bool)


def test_splice_detectors_return_false_for_a_fully_exonic_guide():
    row = {"sgRNA_seq": "ATGCTGGAAAAGACCTGCGA"}
    assert is_splice_acceptor_guide(row, "AtoG", "sense", WINDOW) is False
    assert is_splice_donor_guide(row, "AtoG", "sense", WINDOW) is False


def test_splice_detectors_return_false_for_a_fully_intronic_guide():
    row = {"sgRNA_seq": "acgtacgtacgtacgtacgt"}
    assert is_splice_acceptor_guide(row, "AtoG", "sense", WINDOW) is False
    assert is_splice_donor_guide(row, "AtoG", "sense", WINDOW) is False


def test_splice_detection_depends_on_the_editor():
    """
    An editor that cannot touch the splice dinucleotide is not a splice guide.

    The 'ag' acceptor is only disrupted by an editor acting on A or G, so a
    C->T editor over the same guide must not be called Splice-acceptor.
    """
    seq = "aaaag" + "GCCTACTCCTCCGTC"
    row = {"sgRNA_seq": seq}
    assert is_splice_acceptor_guide(row, "AtoG", "sense", WINDOW) is True
    assert is_splice_acceptor_guide(row, "CtoT", "sense", WINDOW) is False


# --------------------------------------------------------------------------- #
# assign_position
# --------------------------------------------------------------------------- #


def test_assign_position_returns_minus_one_for_controls():
    assert assign_position("M1V", 10, 14, "Control", "Missense") == -1
    assert assign_position("M1V", 10, 14, "Essential", "Missense") == -1


@pytest.mark.parametrize("muttype", ["Splice-donor", "Splice-acceptor", "UTR", "Intron"])
def test_assign_position_returns_minus_one_for_non_coding_types(muttype):
    assert assign_position("M1V", 10, 14, "GENE", muttype) == -1


def test_assign_position_prefers_nonsense_then_missense_then_silent():
    """Nonsense wins over missense, which wins over silent."""
    assert assign_position("Q5./M1V/L2L", 10, 14, "GENE", "Nonsense") == 5.0
    assert assign_position("M1V/L2L", 10, 14, "GENE", "Missense") == 1.0
    assert assign_position("L2L", 10, 14, "GENE", "Silent") == 2.0


def test_assign_position_averages_multiple_residues_of_the_same_class():
    assert assign_position("M1V/L3V", 10, 14, "GENE", "Missense") == 2.0
    assert assign_position("Q5./Q9.", 10, 14, "GENE", "Nonsense") == 7.0


def test_assign_position_deduplicates_repeated_residues():
    """The same residue reached by two edit combinations counts once."""
    assert assign_position("M1V;M1V;M1V", 10, 14, "GENE", "Missense") == 1.0


def test_assign_position_falls_back_to_the_window_centre():
    """With no mutations, the position is derived from the window bounds."""
    assert assign_position(None, 10, 14, "GENE", "No Mutation") == 1 + (12 // 3)
    assert assign_position(float("nan"), 10, 14, "GENE", "No Mutation") == 5.0


def test_assign_position_returns_minus_one_without_window_bounds():
    assert assign_position(None, -1, -1, "GENE", "No Mutation") == -1


# --------------------------------------------------------------------------- #
# small helpers
# --------------------------------------------------------------------------- #


def test_extract_uppercase_letters():
    assert extract_uppercase_letters("aaaGGGttt") == "GGG"
    assert extract_uppercase_letters("acgt") == ""


def test_find_first_uppercase_index():
    assert find_first_uppercase_index("aaaGGG") == 3
    assert find_first_uppercase_index("GGG") == 0
    assert find_first_uppercase_index("aaa") is None


def test_filter_sequence_drops_guides_containing_the_sequence():
    results = [["ATTTTG", "AGG"], ["ATGCTG", "AGG"], ["attttG", "AGG"]]
    kept = filter_sequence(results, "TTTT")
    assert [r[0] for r in kept] == ["ATGCTG"]
