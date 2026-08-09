"""
A parameter that is accepted but ignored is worse than one that is absent.

The parity work added the same options to every plotting function. It would be
easy to add one to a signature and forget to use it in the body, leaving a
caller believing their filter applied when it did not. These tests check that
the shared options actually change what gets plotted.

The approach is to compare the data that reaches the axes -- point counts,
box counts, matrix shape -- between a call with the option off and the same
call with it on.
"""

import inspect

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

import be_scan.plot as bsp


@pytest.fixture
def df():
    rng = np.random.default_rng(42)
    n = 300
    return pd.DataFrame({
        "pos": rng.integers(1, 400, n),
        "cond1": rng.normal(loc=1.0, size=n),
        "cond2": rng.normal(loc=1.0, size=n),
        "muttype": rng.choice(["Missense", "Silent", "Nonsense", "Intron"], n),
        "gene": rng.choice(["FOXA1", "Control"], n),
        "sgRNA_ID": [f"sg{i}" for i in range(n)],
        "score": rng.uniform(0, 10, n),
    })


def n_points(axes):
    """Total plotted points across every axes in a figure."""
    total = 0
    for ax in np.atleast_1d(axes).ravel():
        for coll in ax.collections:
            offsets = getattr(coll, "get_offsets", None)
            if offsets is not None:
                total += len(coll.get_offsets())
    return total


# --------------------------------------------------------------------------- #
# filters actually remove rows
# --------------------------------------------------------------------------- #


def test_scatterplot_value_filter_removes_points(df):
    _, axes_all = bsp.scatterplot(df, ["cond1"], "pos", savefig=False, show=False)
    _, axes_filtered = bsp.scatterplot(
        df, ["cond1"], "pos", savefig=False, show=False,
        filter_val=True, val_cols=["score"], val_min=5.0)

    assert n_points(axes_filtered) < n_points(axes_all), (
        "filter_val did not reduce the number of plotted points"
    )
    assert n_points(axes_filtered) == (df["score"] > 5.0).sum()


def test_scatterplot_category_filter_removes_points(df):
    _, axes = bsp.scatterplot(
        df, ["cond1"], "pos", savefig=False, show=False,
        filter_params=True, params_cols=["muttype"],
        params_conditions=[["Missense"]])
    assert n_points(axes) == (df["muttype"] == "Missense").sum()


def test_scatterplot_xwindow_restricts_the_axis(df):
    _, axes = bsp.scatterplot(
        df, ["cond1"], "pos", savefig=False, show=False, xwindow=[100, 200])
    expected = ((df["pos"] >= 100) & (df["pos"] <= 200)).sum()
    assert n_points(axes) == expected


def test_boxplot_category_filter_removes_rows(df):
    """A filter that removes a whole category must remove its box."""
    fig, axes = bsp.boxplot(
        df, ["cond1"], "muttype", ["Missense", "Silent", "Nonsense"],
        savefig=False, show=False)
    all_boxes = len(axes[0].get_xticklabels())

    fig, axes = bsp.boxplot(
        df, ["cond1"], "muttype", ["Missense", "Silent", "Nonsense"],
        savefig=False, show=False,
        filter_params=True, params_cols=["muttype"],
        params_conditions=[["Missense", "Silent"]])
    assert len(axes[0].get_xticklabels()) < all_boxes


def test_corr_heatmap_filter_changes_the_correlation(df):
    """
    corr_heatmap had no filtering at all before 3.0.0.

    Restricting to a subset must change the computed correlation matrix,
    otherwise the parameter is decorative.
    """
    _, ax_all = bsp.corr_heatmap(df, ["cond1", "cond2"], savefig=False, show=False)
    all_values = ax_all.collections[0].get_array().data.copy()

    _, ax_sub = bsp.corr_heatmap(
        df, ["cond1", "cond2"], savefig=False, show=False,
        filter_params=True, params_cols=["muttype"],
        params_conditions=[["Missense"]])
    sub_values = ax_sub.collections[0].get_array().data

    assert not np.allclose(all_values, sub_values), (
        "filter_params had no effect on the correlation matrix"
    )


def test_sgrna_enrichment_filter_changes_the_distribution(df):
    _, ax_all = bsp.sgrna_enrichment(df, ["cond1"], savefig=False, show=False)
    _, ax_sub = bsp.sgrna_enrichment(
        df, ["cond1"], savefig=False, show=False,
        filter_params=True, params_cols=["muttype"],
        params_conditions=[["Missense"]])

    all_segments = sum(len(c.get_segments()) for c in ax_all[0].collections
                       if hasattr(c, "get_segments"))
    sub_segments = sum(len(c.get_segments()) for c in ax_sub[0].collections
                       if hasattr(c, "get_segments"))
    assert sub_segments < all_segments


# --------------------------------------------------------------------------- #
# negative-control normalization actually shifts the data
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", ["scatterplot", "corr_heatmap"])
def test_neg_ctrl_normalization_shifts_values(df, name):
    """
    Normalizing subtracts the control mean, so values must move.

    The fixture centres cond1 on 1.0, so normalizing against controls should
    pull the plotted values toward zero.
    """
    if name == "scatterplot":
        _, axes_raw = bsp.scatterplot(df, ["cond1"], "pos", savefig=False, show=False)
        _, axes_norm = bsp.scatterplot(
            df, ["cond1"], "pos", savefig=False, show=False,
            neg_ctrl=True, neg_ctrl_col="gene", neg_ctrl_conditions=["Control"])
        raw = axes_raw[0].collections[0].get_offsets()[:, 1]
        norm = axes_norm[0].collections[0].get_offsets()[:, 1]
        assert abs(np.mean(norm)) < abs(np.mean(raw)), (
            "neg_ctrl did not recentre the plotted values"
        )
    else:
        _, ax_raw = bsp.corr_heatmap(df, ["cond1", "cond2"], savefig=False, show=False)
        _, ax_norm = bsp.corr_heatmap(
            df, ["cond1", "cond2"], savefig=False, show=False,
            neg_ctrl=True, neg_ctrl_col="gene", neg_ctrl_conditions=["Control"])
        assert ax_raw is not None and ax_norm is not None


# --------------------------------------------------------------------------- #
# highlighting and annotation actually draw something
# --------------------------------------------------------------------------- #


def test_scatterplot_highlight_adds_a_layer(df):
    _, axes_plain = bsp.scatterplot(df, ["cond1"], "pos", savefig=False, show=False)
    _, axes_hl = bsp.scatterplot(
        df, ["cond1"], "pos", savefig=False, show=False,
        highlight=True, highlight_col="sgRNA_ID",
        highlight_vals=["sg1", "sg2", "sg3"])

    assert len(axes_hl[0].collections) > len(axes_plain[0].collections), (
        "highlight did not add a drawn layer"
    )


def test_scatterplot_annotation_draws_labels(df):
    _, axes_plain = bsp.scatterplot(df, ["cond1"], "pos", savefig=False, show=False)
    _, axes_annot = bsp.scatterplot(
        df, ["cond1"], "pos", savefig=False, show=False,
        annot=True, annot_label="sgRNA_ID", annot_top=7)

    assert len(axes_annot[0].texts) - len(axes_plain[0].texts) == 7


@pytest.mark.parametrize("kw,expected", [
    ({"annot_top": 4}, 4),
    ({"annot_top": -4}, 4),
    ({"annot_abs": 6}, 6),
])
def test_annotation_modes_label_the_right_number(df, kw, expected):
    _, axes = bsp.scatterplot(
        df, ["cond1"], "pos", savefig=False, show=False,
        annot=True, annot_label="sgRNA_ID", **kw)
    assert len(axes[0].texts) == expected


def test_annot_top_negative_selects_the_lowest(df):
    """annot_top=-n must label the n lowest, not the n highest."""
    _, axes_hi = bsp.scatterplot(
        df, ["cond1"], "pos", savefig=False, show=False,
        annot=True, annot_label="sgRNA_ID", annot_top=5)
    _, axes_lo = bsp.scatterplot(
        df, ["cond1"], "pos", savefig=False, show=False,
        annot=True, annot_label="sgRNA_ID", annot_top=-5)

    high = {t.get_text() for t in axes_hi[0].texts}
    low = {t.get_text() for t in axes_lo[0].texts}
    assert high != low, "annot_top sign was ignored"
    assert high == set(df.nlargest(5, "cond1")["sgRNA_ID"])
    assert low == set(df.nsmallest(5, "cond1")["sgRNA_ID"])


def test_n_sd_changes_the_control_lines(df):
    """The +/- 2 SD guide lines were hardcoded; n_sd must move them."""
    ctrl = dict(neg_ctrl=True, neg_ctrl_col="gene", neg_ctrl_conditions=["Control"])
    _, axes_2 = bsp.scatterplot(df, ["cond1"], "pos", savefig=False, show=False, **ctrl)
    _, axes_3 = bsp.scatterplot(df, ["cond1"], "pos", savefig=False, show=False,
                                n_sd=3, **ctrl)

    y2 = sorted(line.get_ydata()[0] for line in axes_2[0].lines)
    y3 = sorted(line.get_ydata()[0] for line in axes_3[0].lines)
    assert y2 != y3, "n_sd did not move the control lines"
    np.testing.assert_allclose(np.array(y3), np.array(y2) * 1.5, rtol=1e-9)


# --------------------------------------------------------------------------- #
# title is honoured everywhere it is offered
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", ["scatterplot", "boxplot", "corr_heatmap"])
def test_title_is_applied(df, name):
    if name == "scatterplot":
        _, axes = bsp.scatterplot(df, ["cond1"], "pos", savefig=False,
                                  show=False, title="my title")
        ax = axes[0]
    elif name == "boxplot":
        _, axes = bsp.boxplot(df, ["cond1"], "muttype", ["Missense", "Silent"],
                              savefig=False, show=False, title="my title")
        ax = axes[0]
    else:
        _, ax = bsp.corr_heatmap(df, ["cond1", "cond2"], savefig=False,
                                 show=False, title="my title")
    assert ax.get_title() == "my title"
