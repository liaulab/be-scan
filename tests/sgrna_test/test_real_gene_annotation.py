"""
End-to-end annotation checks against a real gene.

The synthetic genes in test_annotation_math.py are built to stress specific
junction phases. This module runs the same cross-check against real FOXA1
input, so the suite also covers whatever a genuine UCSC-formatted .fasta
happens to contain: uneven exon sizes, real splice motifs, real UTRs.
"""

import re
import warnings

import pytest

from be_scan.sgrna._genomic_ import DNA_AA_map
from be_scan.sgrna.annotate import annotate
from be_scan.sgrna.generate_library import generate_library

FASTA_NAME = "240401_FOXA1_Input_30.fasta"

EDITORS = [
    pytest.param("A", "G", id="ABE_AtoG"),
    pytest.param("C", "T", id="CBE_CtoT"),
    pytest.param("AC", "GT", id="dual_ACtoGT"),
]


@pytest.fixture(scope="module")
def foxa1_fasta(request):
    path = request.config.rootpath / "tests" / "test_data" / "sgrna" / FASTA_NAME
    if not path.exists():
        pytest.skip(f"{FASTA_NAME} not available")
    return path


@pytest.fixture(scope="module")
def foxa1_library(foxa1_fasta):
    """Generate the guide library once and reuse it across editors."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        guides, gene = generate_library(
            gene_filepath=str(foxa1_fasta),
            gene_name="FOXA1",
            cas_type="SpG",
            edit_from="A",
            edit_to="G",
            window=[4, 8],
            return_df=True,
            save_df=False,
        )
    cds = "".join(gene.exons)
    protein = "".join(DNA_AA_map[cds[i:i + 3]] for i in range(0, len(cds), 3))
    return {"guides": guides, "gene": gene, "cds": cds, "protein": protein}


# --------------------------------------------------------------------------- #
# the gene itself parses correctly
# --------------------------------------------------------------------------- #


def test_foxa1_parses_into_a_clean_reading_frame(foxa1_library):
    cds, protein = foxa1_library["cds"], foxa1_library["protein"]

    assert len(cds) % 3 == 0, "FOXA1 CDS is not a whole number of codons"
    assert cds.startswith("ATG"), "FOXA1 CDS does not start with a start codon"
    assert protein.endswith("."), "FOXA1 CDS does not end with a stop codon"
    assert "." not in protein[:-1], "FOXA1 translation has an internal stop"
    assert len(protein) == 473, f"FOXA1 should be 473 residues, got {len(protein)}"


def test_foxa1_metadata_is_extracted(foxa1_library):
    gene = foxa1_library["gene"]
    assert gene.chrID == "chr14"
    assert gene.strand == "minus"
    assert gene.intron_len == 30
    assert len(gene.exons) == 2


# --------------------------------------------------------------------------- #
# every reported mutation agrees with the translated CDS
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("edit_from,edit_to", EDITORS)
def test_every_mutation_call_matches_the_translated_protein(
    foxa1_library, edit_from, edit_to
):
    """
    For every 'X###Y' emitted, residue ### of the translated CDS must be X.

    This catches frame errors and off-by-one residue numbering across the
    whole library at once.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annotated = annotate(
            guides_file=foxa1_library["guides"],
            edit_from=edit_from,
            edit_to=edit_to,
            exons=foxa1_library["gene"].exons,
            window=[4, 8],
            return_df=True,
            save_df=False,
        )

    protein = foxa1_library["protein"]
    pre = f"{edit_from}to{edit_to}"

    checked, mismatches = 0, []
    for _, row in annotated.iterrows():
        cell = row[f"{pre}_mutations"]
        if not isinstance(cell, str) or not cell:
            continue
        for mut in re.split(r"[;/]", cell):
            match = re.fullmatch(r"([A-Z.])(\d+)([A-Z.])", mut)
            if not match:
                continue
            wt, position = match.group(1), int(match.group(2))
            checked += 1
            if not 1 <= position <= len(protein):
                mismatches.append((row["sgRNA_ID"], mut, "position out of range"))
            elif protein[position - 1] != wt:
                mismatches.append(
                    (row["sgRNA_ID"], mut, f"protein[{position}] is {protein[position - 1]}")
                )

    assert checked > 500, f"only {checked} mutation calls checked"
    assert not mismatches, (
        f"{len(mismatches)}/{checked} mutation calls disagree with the "
        f"translated CDS:\n" + "\n".join(f"  {m}" for m in mismatches[:15])
    )


# --------------------------------------------------------------------------- #
# annotation output is well formed
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("edit_from,edit_to", EDITORS)
def test_annotation_adds_the_expected_columns(foxa1_library, edit_from, edit_to):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annotated = annotate(
            guides_file=foxa1_library["guides"],
            edit_from=edit_from, edit_to=edit_to,
            exons=foxa1_library["gene"].exons,
            window=[4, 8], return_df=True, save_df=False,
        )

    pre = f"{edit_from}to{edit_to}"
    for suffix in ("win_overlap", "mutations", "muttypes", "muttype", "pos"):
        assert f"{pre}_{suffix}" in annotated.columns

    # intermediate columns must be dropped again
    for suffix in ("target_CDS", "codon_window", "residue_window", "target_windowpos"):
        assert f"{pre}_{suffix}" not in annotated.columns

    # the input columns must survive untouched
    for col in ("sgRNA_seq", "starting_frame", "sgRNA_strand", "gene_pos", "exon"):
        assert col in annotated.columns


def test_muttype_values_come_from_the_documented_set(foxa1_library):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annotated = annotate(
            guides_file=foxa1_library["guides"],
            edit_from="A", edit_to="G",
            exons=foxa1_library["gene"].exons,
            window=[4, 8], return_df=True, save_df=False,
        )

    allowed = {
        "Nonsense", "Missense", "Silent", "Splice-acceptor", "Splice-donor",
        "Intron", "UTR", "No Mutation", "Flank",
    }
    seen = set(annotated["AtoG_muttype"].dropna().unique())
    assert seen <= allowed, f"unexpected muttype values: {seen - allowed}"


def test_win_overlap_agrees_with_the_guide_sequence(foxa1_library):
    """
    Exon / Intron / Exon-Intron must match the case of the window bases.

    Sense guides read coding_seq[3:8]; antisense guides read the mirrored
    window, so the two are checked separately.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annotated = annotate(
            guides_file=foxa1_library["guides"],
            edit_from="A", edit_to="G",
            exons=foxa1_library["gene"].exons,
            window=[4, 8], return_df=True, save_df=False,
        )

    for _, row in annotated.iterrows():
        coding = row["coding_seq"]
        if row["sgRNA_strand"] == "sense":
            window_seq = coding[3:8]
        else:
            window_seq = coding[20 - 8:20 - 4 + 1]

        if window_seq.isupper():
            expected = "Exon"
        elif window_seq.islower():
            expected = "Intron"
        else:
            expected = "Exon/Intron"

        assert row["AtoG_win_overlap"] == expected, (
            f"{row['sgRNA_ID']} ({row['sgRNA_strand']}): window '{window_seq}' "
            f"should be {expected}, got {row['AtoG_win_overlap']}"
        )


def test_edit_base_count_matches_the_window(foxa1_library):
    """The X_count column counts editable bases in window positions 4-8."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annotated = annotate(
            guides_file=foxa1_library["guides"],
            edit_from="A", edit_to="G",
            exons=foxa1_library["gene"].exons,
            window=[4, 8], return_df=True, save_df=False,
        )

    for _, row in annotated.iterrows():
        expected = row["sgRNA_seq"][3:8].count("A")
        assert row["A_count"] == expected, (
            f"{row['sgRNA_ID']}: A_count is {row['A_count']}, "
            f"window '{row['sgRNA_seq'][3:8]}' has {expected}"
        )


def test_positions_are_within_the_protein_or_flagged(foxa1_library):
    """AtoG_pos is either -1 (not applicable) or a residue in range."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        annotated = annotate(
            guides_file=foxa1_library["guides"],
            edit_from="A", edit_to="G",
            exons=foxa1_library["gene"].exons,
            window=[4, 8], return_df=True, save_df=False,
        )

    n_protein = len(foxa1_library["protein"])
    positions = annotated["AtoG_pos"].dropna()
    out_of_range = positions[(positions != -1) & ((positions < 1) | (positions > n_protein))]
    assert out_of_range.empty, f"positions outside 1..{n_protein}: {list(out_of_range)[:10]}"


# --------------------------------------------------------------------------- #
# determinism
# --------------------------------------------------------------------------- #


def test_annotation_is_deterministic(foxa1_library):
    """
    Two runs over the same input must agree.

    annotate joins mutation strings out of a set, so ordering within a cell is
    not stable; the comparison is made order-insensitive before checking.
    """
    def run():
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return annotate(
                guides_file=foxa1_library["guides"],
                edit_from="A", edit_to="G",
                exons=foxa1_library["gene"].exons,
                window=[4, 8], return_df=True, save_df=False,
            )

    first, second = run(), run()
    assert list(first["sgRNA_ID"]) == list(second["sgRNA_ID"])

    def normalize(cell):
        if not isinstance(cell, str):
            return cell
        return frozenset(frozenset(c.split("/")) for c in cell.split(";"))

    a = first["AtoG_mutations"].map(normalize)
    b = second["AtoG_mutations"].map(normalize)
    assert a.equals(b), "annotation is not reproducible across runs"

    for col in ("AtoG_muttype", "AtoG_pos", "AtoG_win_overlap"):
        assert first[col].equals(second[col]), f"{col} differs between runs"
