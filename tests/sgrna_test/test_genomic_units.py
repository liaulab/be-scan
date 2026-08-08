"""Unit tests for be_scan.sgrna._genomic_ -- the primitives everything rests on."""

import re

import pytest

from be_scan.sgrna._genomic_ import (
    DNA_AA_map,
    bases,
    cas_key,
    complement,
    complements,
    DNA_to_AA,
    process_PAM,
    protein_to_AAseq,
    rev_complement,
    first_capital_index,
    last_capital_index,
)


# --------------------------------------------------------------------------- #
# the codon table itself
# --------------------------------------------------------------------------- #


def test_codon_table_is_complete():
    """All 64 codons present, exactly 3 stops, exactly 1 methionine codon."""
    assert len(DNA_AA_map) == 64
    expected = {a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT"}
    assert set(DNA_AA_map) == expected
    stops = [c for c, aa in DNA_AA_map.items() if aa == "."]
    assert sorted(stops) == ["TAA", "TAG", "TGA"]
    assert [c for c, aa in DNA_AA_map.items() if aa == "M"] == ["ATG"]


def test_codon_table_matches_standard_genetic_code():
    """Spot-check against the standard code, including the 6-fold families."""
    for codon in ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"):
        assert DNA_AA_map[codon] == "L"
    for codon in ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"):
        assert DNA_AA_map[codon] == "S"
    for codon in ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"):
        assert DNA_AA_map[codon] == "R"
    assert DNA_AA_map["TGG"] == "W"  # only Trp codon
    assert DNA_AA_map["ATG"] == "M"


# --------------------------------------------------------------------------- #
# DNA_to_AA
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "seq,expected",
    [
        ("ATG", "M"),
        ("ATGTTAGGA", "MLG"),
        ("TAA", "."),
        ("TAG", "."),
        ("TGA", "."),
        ("", ""),
    ],
)
def test_dna_to_aa_uppercase(seq, expected):
    assert DNA_to_AA(seq) == expected


def test_dna_to_aa_lowercase_is_translated_when_upper_is_true():
    """upper=True uppercases first, so intronic sequence translates normally."""
    assert DNA_to_AA("atgtta", upper=True) == "ML"


def test_dna_to_aa_lowercase_becomes_underscore_when_upper_is_false():
    """
    upper=False is how annotation marks non-coding codons.

    A codon is only translated if it is entirely uppercase; anything with a
    lowercase base becomes '_', which format_mutation then skips.
    """
    assert DNA_to_AA("atgTTA", upper=False) == "_L"
    assert DNA_to_AA("ATGtta", upper=False) == "M_"
    assert DNA_to_AA("atgtta", upper=False) == "__"
    # a codon that is only partly lowercase is still non-coding
    assert DNA_to_AA("aTGTTA", upper=False) == "_L"


@pytest.mark.parametrize("bad", ["AT", "ATGT", "ATGTT"])
def test_dna_to_aa_rejects_non_multiple_of_three(bad):
    with pytest.raises(AssertionError):
        DNA_to_AA(bad)


@pytest.mark.parametrize("bad", ["ATN", "AT-", "AT ", "ATX"])
def test_dna_to_aa_rejects_non_nucleotides(bad):
    with pytest.raises(AssertionError):
        DNA_to_AA(bad)


def test_dna_to_aa_rejects_non_string():
    with pytest.raises(AssertionError):
        DNA_to_AA(["A", "T", "G"])


# --------------------------------------------------------------------------- #
# complement / rev_complement
# --------------------------------------------------------------------------- #


def test_complement_pairs():
    assert complements == {
        "A": "T", "T": "A", "G": "C", "C": "G",
        "a": "t", "t": "a", "g": "c", "c": "g",
    }


@pytest.mark.parametrize(
    "seq,expected",
    [
        ("ATGC", "GCAT"),
        ("AAAA", "TTTT"),
        ("acgt", "acgt"),
        ("ATGCatgc", "gcatGCAT"),
        ("", ""),
    ],
)
def test_rev_complement(seq, expected):
    assert rev_complement(complements, seq) == expected


@pytest.mark.parametrize("seq", ["ATGC", "acgtACGT", "GGGCCC", "A"])
def test_rev_complement_is_an_involution(seq):
    assert rev_complement(complements, rev_complement(complements, seq)) == seq


def test_rev_complement_preserves_case():
    """Case carries the exon/intron distinction and must survive the flip."""
    out = rev_complement(complements, "ATGCatgc")
    assert out == "gcatGCAT"
    assert sum(c.isupper() for c in out) == 4


def test_complement_is_not_reversed():
    assert complement(complements, "ATGC") == "TACG"
    assert complement(complements, "ATGC") == rev_complement(complements, "ATGC")[::-1]


@pytest.mark.parametrize("bad", ["ATN", "xyz", "AT GC"])
def test_complement_functions_reject_bad_input(bad):
    with pytest.raises(AssertionError):
        rev_complement(complements, bad)
    with pytest.raises(AssertionError):
        complement(complements, bad)


# --------------------------------------------------------------------------- #
# process_PAM
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize(
    "pam,matching,nonmatching",
    [
        ("NGG", ["AGG", "TGG", "CGG", "GGG"], ["AGA", "AAG", "TTT"]),
        ("NGN", ["AGA", "TGC", "CGG", "GGT"], ["AAA", "ATC"]),
        ("NNN", ["AAA", "TGC", "GGG"], []),
        ("NGA", ["AGA", "TGA"], ["AGG", "AAA"]),
        ("NRN", ["AAA", "AGT", "TAC"], ["ACA", "ATG"]),
        ("NYN", ["ACA", "ATG", "TCC"], ["AAA", "AGT"]),
    ],
)
def test_process_pam_matches_expected_sequences(pam, matching, nonmatching):
    regex = process_PAM(pam)
    for seq in matching:
        assert regex.match(seq), f"{pam} should match {seq}"
    for seq in nonmatching:
        assert not regex.match(seq), f"{pam} should not match {seq}"


def test_process_pam_is_case_insensitive_on_input_and_target():
    assert process_PAM("ngg").match("AGG")
    assert process_PAM("NGG").match("agg")


@pytest.mark.parametrize("bad", ["NGX", "123", "NG-"])
def test_process_pam_rejects_invalid_characters(bad):
    with pytest.raises(AssertionError):
        process_PAM(bad)


def test_every_named_cas_has_a_processable_pam():
    for name, pam in cas_key.items():
        regex = process_PAM(pam)
        assert isinstance(regex, re.Pattern), f"{name} produced no regex"
        assert len(pam) == 3, f"{name} PAM {pam} is not 3 nt"


def test_bases_constant():
    assert bases == "ACGT"


# --------------------------------------------------------------------------- #
# protein_to_AAseq
# --------------------------------------------------------------------------- #


def test_protein_to_aaseq_is_one_indexed_with_trailing_stop(tmp_path):
    """
    The dict is 1-indexed and gains a synthetic stop one past the end.

    annotate uses it as amino_acid_seq[start+i], so the indexing convention
    matters: residue 1 must be the first amino acid.
    """
    fasta = tmp_path / "p.fasta"
    fasta.write_text(">sp|TEST|test protein\nMLGTV\nKMEGH\n")

    seq = protein_to_AAseq(str(fasta))
    assert seq[1] == "M"
    assert seq[2] == "L"
    assert seq[10] == "H"
    assert seq[11] == "."  # synthetic stop past the last residue
    assert 0 not in seq


# --------------------------------------------------------------------------- #
# capital index helpers
# --------------------------------------------------------------------------- #


def test_first_and_last_capital_index():
    assert first_capital_index("aaaATGttt") == 3
    assert last_capital_index("aaaATGttt") == 5
    assert first_capital_index("ATG") == 0
    assert last_capital_index("ATG") == 2


@pytest.mark.parametrize("fn", [first_capital_index, last_capital_index])
def test_capital_index_raises_when_all_lowercase(fn):
    with pytest.raises(Exception):
        fn("acgtacgt")
