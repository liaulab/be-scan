"""
The sgRNA spacer length is a parameter, not a constant.

be_scan assumed 20 nt spacers in six places: GeneForCRISPR's guide_len
attribute, the UTR bounds in find_all_guides, the reverse-guide PAM split,
generate_library's window assertion, annotate_intron_exon_anti's mirrored
slice, and count_reads' read filter. Cas12a and truncated-spacer libraries
are not 20 nt, so all of them now take guide_len, defaulting to 20.

These tests check two things: that guide_len=20 is byte-identical to the old
behavior, and that other lengths actually produce guides of that length.
"""

import warnings

import pytest

from be_scan.sgrna._gene_ import GeneForCRISPR
from be_scan.sgrna.generate_library import generate_library

from tests.sgrna_test import synthetic_gene as syn


@pytest.fixture
def gene_file(tmp_path):
    return syn.write_gene(tmp_path, syn.ALL_PHASE_LENGTHS)


def build(gene_file, guide_len=20, cas_type="SpG", window=(4, 8), **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return generate_library(
            gene_filepath=str(gene_file["path"]),
            gene_name="SYN", cas_type=cas_type,
            edit_from="A", edit_to="G", window=list(window),
            return_df=True, save_df=False, guide_len=guide_len, **kw,
        )


# --------------------------------------------------------------------------- #
# the default is unchanged
# --------------------------------------------------------------------------- #


def test_guide_len_defaults_to_twenty():
    import inspect
    assert inspect.signature(generate_library).parameters["guide_len"].default == 20
    assert inspect.signature(GeneForCRISPR.__init__).parameters["guide_len"].default == 20

    from be_scan.analysis.count_reads import count_reads
    from be_scan.analysis.batch_process import batch_process
    assert inspect.signature(count_reads).parameters["guide_len"].default == 20
    assert inspect.signature(batch_process).parameters["guide_len"].default == 20


def test_explicit_twenty_matches_the_default(gene_file):
    """Passing guide_len=20 must be indistinguishable from not passing it."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        implicit, _ = generate_library(
            gene_filepath=str(gene_file["path"]), gene_name="SYN", cas_type="SpG",
            edit_from="A", edit_to="G", window=[4, 8],
            return_df=True, save_df=False)
    explicit, _ = build(gene_file, guide_len=20)

    assert implicit.shape == explicit.shape
    assert list(implicit.columns) == list(explicit.columns)
    for col in implicit.columns:
        assert implicit[col].equals(explicit[col]), f"{col} differs"


# --------------------------------------------------------------------------- #
# other lengths actually work
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("guide_len", [18, 19, 20, 21, 22, 23])
def test_generated_guides_have_the_requested_length(gene_file, guide_len):
    guides, gene = build(gene_file, guide_len=guide_len)

    assert len(guides) > 0, f"no guides generated at guide_len={guide_len}"
    lengths = guides["sgRNA_seq"].str.len().unique()
    assert set(lengths) == {guide_len}, (
        f"guide_len={guide_len} produced spacers of length {sorted(lengths)}"
    )
    assert gene.guide_len == guide_len


@pytest.mark.parametrize("guide_len", [18, 20, 22])
def test_pam_length_is_consistent(gene_file, guide_len):
    """PAM_seq holds whatever the n-mer has left after the spacer."""
    guides, _ = build(gene_file, guide_len=guide_len, cas_type="SpG")
    pam_lengths = guides["PAM_seq"].str.len().unique()
    assert set(pam_lengths) == {3}, f"PAM lengths {sorted(pam_lengths)}"


@pytest.mark.parametrize("pam", ["NGG", "NGN", "NNN"])
def test_custom_pam_length_is_respected(gene_file, pam):
    """The searched n-mer is guide_len + len(PAM), not a hardcoded 23."""
    guides, _ = build(gene_file, guide_len=20, PAM=pam)
    assert set(guides["sgRNA_seq"].str.len().unique()) == {20}
    assert set(guides["PAM_seq"].str.len().unique()) == {len(pam)}


def test_both_strands_respect_the_guide_length(gene_file):
    """The reverse-guide split used a hardcoded PAM length of 3."""
    guides, _ = build(gene_file, guide_len=22)
    for strand in ("sense", "antisense"):
        subset = guides[guides["sgRNA_strand"] == strand]
        assert len(subset) > 0, f"no {strand} guides"
        assert set(subset["sgRNA_seq"].str.len().unique()) == {22}, (
            f"{strand} guides are not 22 nt"
        )


# --------------------------------------------------------------------------- #
# the window assertion follows guide_len
# --------------------------------------------------------------------------- #


def test_window_may_span_the_whole_guide(gene_file):
    """A window up to guide_len is accepted, not capped at 20."""
    guides, _ = build(gene_file, guide_len=23, window=(4, 22))
    assert len(guides) > 0


def test_window_beyond_the_guide_is_rejected(gene_file):
    with pytest.raises(AssertionError, match="fit within"):
        build(gene_file, guide_len=18, window=(4, 19))


def test_window_beyond_a_20mer_is_still_rejected(gene_file):
    """The previous hardcoded limit is preserved at the default length."""
    with pytest.raises(AssertionError, match="fit within"):
        build(gene_file, guide_len=20, window=(4, 21))


# --------------------------------------------------------------------------- #
# count_reads
# --------------------------------------------------------------------------- #


def _fastq(path, reads, key="CGAAACACC", key_rev="GTTTGAGA"):
    """Minimal FASTQ whose reads carry a spacer between the two keys."""
    lines = []
    for i, guide in enumerate(reads):
        seq = ("N" * 10) + key + guide + key_rev + ("N" * 10)
        lines += [f"@read{i}", seq, "+", "I" * len(seq)]
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.mark.parametrize("guide_len", [18, 20, 22])
def test_count_reads_counts_guides_of_the_given_length(tmp_path, guide_len):
    import pandas as pd
    from be_scan.analysis.count_reads import count_reads

    spacers = ["ACGT" * 6, "TGCA" * 6, "GGCC" * 6]
    spacers = [s[:guide_len] for s in spacers]

    _fastq(tmp_path / "r.fastq", spacers * 4)
    lib = pd.DataFrame({"sgRNA_seq": spacers})
    sheet = pd.DataFrame({
        "fastq_file": ["r.fastq"], "counts_file": ["c.csv"],
        "noncounts_file": ["n.csv"], "stats_file": ["s.txt"],
        "condition": ["t0"],
    })

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = count_reads(
            sample_sheet=sheet, annotated_lib=lib,
            in_dir=str(tmp_path), out_dir=str(tmp_path),
            return_df=True, save_files=False, guide_len=guide_len,
        )

    assert out["t0"].sum() == 12, (
        f"guide_len={guide_len}: counted {out['t0'].sum()} of 12 reads"
    )


def test_count_reads_rejects_the_wrong_length(tmp_path):
    """A library of 22-mers read with the default guide_len=20 counts nothing."""
    import pandas as pd
    from be_scan.analysis.count_reads import count_reads

    spacers = [("ACGT" * 6)[:22], ("TGCA" * 6)[:22]]
    _fastq(tmp_path / "r.fastq", spacers * 3)
    lib = pd.DataFrame({"sgRNA_seq": spacers})
    sheet = pd.DataFrame({
        "fastq_file": ["r.fastq"], "counts_file": ["c.csv"],
        "noncounts_file": ["n.csv"], "stats_file": ["s.txt"],
        "condition": ["t0"],
    })

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        default = count_reads(
            sample_sheet=sheet, annotated_lib=lib, in_dir=str(tmp_path),
            out_dir=str(tmp_path), return_df=True, save_files=False)
        matched = count_reads(
            sample_sheet=sheet, annotated_lib=lib, in_dir=str(tmp_path),
            out_dir=str(tmp_path), return_df=True, save_files=False, guide_len=22)

    assert default["t0"].sum() == 0, "22-mers should not match at guide_len=20"
    assert matched["t0"].sum() == 6, "22-mers should match at guide_len=22"
