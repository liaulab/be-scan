"""
Regression tests for defects fixed in 3.0.0.

Each of these was reachable on a default or documented code path. They are
grouped here rather than scattered so the fixes stay fixed.
"""

import warnings

import numpy as np
import pandas as pd
import pytest


# --------------------------------------------------------------------------- #
# B1: calc_controls raised UnboundLocalError on every call
# --------------------------------------------------------------------------- #


def _conditions_frame(n=30, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "sgRNA_seq": [f"g{i}" for i in range(n)],
        "gene": (["FOXA1"] * (n - 10)) + (["Control"] * 10),
        "cond1_LFCminusControl_avg": rng.normal(size=n),
        "cond2_LFCminusControl_avg": rng.normal(size=n),
    })


def test_calc_controls_runs_with_a_dataframe(tmp_path):
    """
    This tested `isinstance(df_conds, ...)` before df_conds existed.

    Every call raised UnboundLocalError, which also made the stats branch of
    batch_process dead code.
    """
    from be_scan.analysis.calc_controls import calc_controls

    calc_controls(
        conditions=_conditions_frame(),
        stats_comparisons=["cond1", "cond2"],
        neg_ctrl_col="gene", neg_ctrl_conditions=["Control"],
        out_dir=str(tmp_path), out_file="stats.txt", save=True,
    )
    text = (tmp_path / "stats.txt").read_text()
    assert "For comparison cond1_LFCminusControl_avg" in text
    assert "Standard deviation is" in text


def test_calc_controls_runs_with_a_csv_path(tmp_path):
    from be_scan.analysis.calc_controls import calc_controls

    csv = tmp_path / "conditions.csv"
    _conditions_frame().to_csv(csv, index=False)

    calc_controls(
        conditions="conditions.csv",
        stats_comparisons=["cond1"],
        neg_ctrl_col="gene", neg_ctrl_conditions=["Control"],
        in_dir=str(tmp_path), out_dir=str(tmp_path),
        out_file="stats.txt", save=True,
    )
    assert (tmp_path / "stats.txt").exists()


def test_calc_controls_does_not_mutate_the_caller_frame():
    from be_scan.analysis.calc_controls import calc_controls

    df = _conditions_frame()
    before = df.copy()
    calc_controls(
        conditions=df, stats_comparisons=["cond1"],
        neg_ctrl_col="gene", neg_ctrl_conditions=["Control"],
        save=False,
    )
    pd.testing.assert_frame_equal(df, before)


# --------------------------------------------------------------------------- #
# B7: count_reads named the wrong file in its sample-sheet error
# --------------------------------------------------------------------------- #


def test_count_reads_error_names_the_sample_sheet(tmp_path):
    from be_scan.analysis.count_reads import count_reads

    bad_sheet = pd.DataFrame({"fastq_file": ["x.fastq"]})  # missing columns
    lib = pd.DataFrame({"sgRNA_seq": ["ACGT" * 5]})

    with pytest.raises(Exception) as excinfo:
        count_reads(sample_sheet=bad_sheet, annotated_lib=lib,
                    in_dir=str(tmp_path), out_dir=str(tmp_path))

    message = str(excinfo.value)
    assert "missing column" in message
    assert "sgRNA_seq" not in message, (
        "the message pointed at annotated_lib while validating the sample sheet"
    )


# --------------------------------------------------------------------------- #
# B8: annotate ignored its own column-name parameters
# --------------------------------------------------------------------------- #


def test_annotate_honours_a_custom_sequence_column(tmp_path):
    """
    seq_col is documented, but three call sites used the literal 'sgRNA_seq'.

    A caller whose frame lacks that literal column got a KeyError; one whose
    frame had both got the wrong column silently.
    """
    from be_scan.sgrna.annotate import annotate
    from be_scan.sgrna.generate_library import generate_library
    from tests.sgrna_test import synthetic_gene as syn

    info = syn.write_gene(tmp_path, syn.ALL_PHASE_LENGTHS)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        guides, gene = generate_library(
            gene_filepath=str(info["path"]), gene_name="SYN", cas_type="SpG",
            edit_from="A", edit_to="G", window=[4, 8],
            return_df=True, save_df=False)

        # rename every column the caller is allowed to rename
        renamed = guides.rename(columns={
            "sgRNA_seq": "spacer",
            "starting_frame": "frame",
            "sgRNA_strand": "strand",
            "gene_pos": "cds_pos",
            "windowstart_pos": "win_start",
            "windowend_pos": "win_end",
            "gene": "gene_name",
        })
        assert "sgRNA_seq" not in renamed.columns

        custom = annotate(
            guides_file=renamed, edit_from="A", edit_to="G",
            exons=gene.exons, window=[4, 8],
            seq_col="spacer", frame_col="frame", strand_col="strand",
            gene_pos_col="cds_pos", window_start_col="win_start",
            window_end_col="win_end", gene_col="gene_name",
            return_df=True, save_df=False)

        default = annotate(
            guides_file=guides, edit_from="A", edit_to="G",
            exons=gene.exons, window=[4, 8],
            return_df=True, save_df=False)

    assert len(custom) == len(default)
    for col in ("AtoG_mutations", "AtoG_muttype", "AtoG_pos", "AtoG_win_overlap"):
        np.testing.assert_array_equal(
            custom[col].fillna("").to_numpy(),
            default[col].fillna("").to_numpy(),
            err_msg=f"{col} differs between custom and default column names")


def test_annotate_gene_col_defaults_to_gene():
    import inspect
    from be_scan.sgrna.annotate import annotate
    assert inspect.signature(annotate).parameters["gene_col"].default == "gene"


# --------------------------------------------------------------------------- #
# B3: DomainOpts.dom_setting was iterated as a list of spans
# --------------------------------------------------------------------------- #


def test_resolve_domains_unwraps_a_gene_keyed_mapping():
    """
    Iterating the mapping yielded its keys, so dom["start"] raised TypeError.

    resolve_domains handles the documented gene-keyed form, a bare list, and
    the empty default.
    """
    from be_scan.figure_plot.figure_classes import DomainOpts, resolve_domains

    spans = [{"start": 10, "end": 50}, {"start": 80, "end": 120}]

    # gene-keyed, selected explicitly
    opts = DomainOpts(dom_setting={"FOXA1": spans, "EZH2": []}, gene="FOXA1")
    assert resolve_domains(opts) == spans

    # gene-keyed with a single gene and no explicit selection
    assert resolve_domains(DomainOpts(dom_setting={"FOXA1": spans})) == spans

    # a bare list is accepted too
    assert resolve_domains(DomainOpts(dom_setting=spans)) == spans

    # the default is empty, not an error
    assert resolve_domains(DomainOpts()) == []


def test_resolve_domains_picks_the_gene_present_in_the_data():
    from be_scan.figure_plot.figure_classes import DomainOpts, resolve_domains

    spans = [{"start": 1, "end": 9}]
    opts = DomainOpts(dom_setting={"FOXA1": spans, "EZH2": [{"start": 2, "end": 3}]})
    df = pd.DataFrame({"Gene Symbol": ["FOXA1"] * 5})
    assert resolve_domains(opts, df) == spans


def test_scatterplot_figure_draws_domains(tmp_path):
    """The end-to-end path that used to raise TypeError."""
    import matplotlib
    matplotlib.use("Agg")
    from be_scan.figure_plot import scatterplot_figure
    from be_scan.figure_plot.figure_classes import DomainOpts, OutputOpts

    df = pd.DataFrame({
        "pos": np.arange(1, 101),
        "score": np.random.default_rng(0).normal(size=100),
        "Gene Symbol": ["FOXA1"] * 100,
    })
    fig, axes = scatterplot_figure(
        df, ["score"], "pos",
        domain=DomainOpts(dom_setting={"FOXA1": [{"start": 20, "end": 40}]},
                          gene="FOXA1"),
        output=OutputOpts(show=False, save=False),
    )
    assert fig is not None


# --------------------------------------------------------------------------- #
# B4: plot.scatterplot called IPython's display() unconditionally
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("annot_kw", [{"annot_top": 5}, {"annot_abs": 5}])
def test_scatterplot_annotation_works_outside_a_notebook(tmp_path, annot_kw):
    """
    A bare `display(...)` only resolves inside IPython.

    Any script using annot_top or annot_abs raised NameError right after
    drawing the labels.
    """
    import matplotlib
    matplotlib.use("Agg")
    from be_scan.plot import scatterplot

    rng = np.random.default_rng(0)
    df = pd.DataFrame({
        "pos": np.arange(1, 61),
        "cond": rng.normal(size=60),
        "label": [f"M{i}V" for i in range(1, 61)],
    })

    scatterplot(df_filepath=df, comparisons=["cond"], x_column="pos",
                annot=True, annot_label="label",
                savefig=False, show=False, out_dir=str(tmp_path), **annot_kw)


# --------------------------------------------------------------------------- #
# B2: correlation_scatterplot_figure raised NameError when not saving
# --------------------------------------------------------------------------- #


def test_correlation_scatterplot_runs_without_saving(tmp_path):
    """plt.close(fig_legend) sat outside the block that created fig_legend."""
    import matplotlib
    matplotlib.use("Agg")
    from be_scan.figure_plot import correlation_scatterplot_figure
    from be_scan.figure_plot.figure_classes import OutputOpts

    rng = np.random.default_rng(0)
    df = pd.DataFrame({
        "x": rng.normal(size=80),
        "y": rng.normal(size=80),
        "mut_type": rng.choice(["Missense", "Silent", "Splice"], size=80),
    })

    correlation_scatterplot_figure(
        df, "x", "y", legend_order=["Missense", "Silent", "Splice"],
        hue_col="mut_type",
        palette={"Missense": "#fb8072", "Silent": "#80b1d3", "Splice": "#fdb462"},
        output=OutputOpts(save=False, show=False),
    )
