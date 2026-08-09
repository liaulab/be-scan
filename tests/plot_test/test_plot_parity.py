"""
Every be_scan.plot function should offer the same capabilities.

Before 3.0.0 the plotting functions had drifted: scatterplot could filter,
colour, annotate and normalize to negative controls; corr_heatmap could do
none of those; sgrna_enrichment was the only one that could highlight a
subset. Two functions taking the same parameter did not always give it the
same meaning.

These tests pin the shared surface down so it cannot drift again. They are
deliberately signature-driven: a new plotting function that forgets the common
parameters fails here rather than at review time.
"""

import inspect

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

import be_scan.plot as bsp


# --------------------------------------------------------------------------- #
# what "parity" means
# --------------------------------------------------------------------------- #

STATIC_PLOTS = [
    "scatterplot", "boxplot", "corr_jointplot", "corr_heatmap",
    "sgrna_enrichment",
]

INTERACTIVE_PLOTS = [
    "interactive_scatter", "interactive_boxplot", "interactive_corr_jointplot",
    "interactive_corr_heatmap",
]

# capability -> the parameters that implement it
CAPABILITIES = {
    "data input":       ["df_filepath"],
    "negative control": ["neg_ctrl", "neg_ctrl_col", "neg_ctrl_conditions"],
    "value filter":     ["filter_val", "val_cols", "val_min"],
    "category filter":  ["filter_params", "params_cols", "params_conditions"],
    "output":           ["savefig", "show", "out_name", "out_dir"],
    "labels":           ["title"],
}

STATIC_ONLY = {
    "figure resolution": ["dpi"],
    "output format":     ["out_type"],
}


def params(fn):
    return set(inspect.signature(fn).parameters)


@pytest.mark.parametrize("name", STATIC_PLOTS)
@pytest.mark.parametrize("capability", sorted(CAPABILITIES))
def test_static_plots_share_every_capability(name, capability):
    fn = getattr(bsp, name)
    missing = set(CAPABILITIES[capability]) - params(fn)
    assert not missing, (
        f"be_scan.plot.{name} is missing {sorted(missing)}, so it cannot do "
        f"'{capability}' while its sibling plots can"
    )


@pytest.mark.parametrize("name", STATIC_PLOTS)
@pytest.mark.parametrize("capability", sorted(STATIC_ONLY))
def test_static_plots_share_output_controls(name, capability):
    fn = getattr(bsp, name)
    missing = set(STATIC_ONLY[capability]) - params(fn)
    assert not missing, f"be_scan.plot.{name} is missing {sorted(missing)}"


@pytest.mark.parametrize("name", INTERACTIVE_PLOTS)
def test_interactive_plots_share_the_core(name):
    fn = getattr(bsp, name)
    required = {"df_filepath", "savefig", "show", "out_name", "out_dir", "title"}
    missing = required - params(fn)
    assert not missing, f"be_scan.plot.{name} is missing {sorted(missing)}"


@pytest.mark.parametrize("name", ["scatterplot", "boxplot", "corr_jointplot",
                                  "sgrna_enrichment"])
def test_plots_that_show_points_can_highlight(name):
    """Highlighting a subset of guides was only available on sgrna_enrichment."""
    fn = getattr(bsp, name)
    missing = {"highlight", "highlight_col", "highlight_vals"} - params(fn)
    assert not missing, f"be_scan.plot.{name} is missing {sorted(missing)}"


@pytest.mark.parametrize("name", ["scatterplot", "boxplot", "corr_jointplot",
                                  "sgrna_enrichment"])
def test_plots_that_show_points_can_colour_by_category(name):
    fn = getattr(bsp, name)
    missing = {"include_hue", "hue_col", "pal"} - params(fn)
    assert not missing, f"be_scan.plot.{name} is missing {sorted(missing)}"


@pytest.mark.parametrize("name", ["scatterplot", "boxplot"])
def test_plots_with_comparisons_can_annotate(name):
    fn = getattr(bsp, name)
    missing = {"annot", "annot_label", "annot_top", "annot_cutoff",
               "annot_abs"} - params(fn)
    assert not missing, f"be_scan.plot.{name} is missing {sorted(missing)}"


# --------------------------------------------------------------------------- #
# the shared parameters actually behave the same way
# --------------------------------------------------------------------------- #


@pytest.fixture
def df():
    rng = np.random.default_rng(0)
    n = 200
    return pd.DataFrame({
        "pos": rng.integers(1, 400, n),
        "cond1": rng.normal(size=n),
        "cond2": rng.normal(size=n),
        "muttype": rng.choice(["Missense", "Silent", "Nonsense"], n),
        "gene": rng.choice(["FOXA1", "Control"], n),
        "sgRNA_ID": [f"sg{i}" for i in range(n)],
        "score": rng.uniform(0, 10, n),
    })


def call(name, df, **kw):
    """Invoke a plot with its required positional arguments plus kw."""
    base = dict(savefig=False, show=False)
    base.update(kw)
    fn = getattr(bsp, name)
    if name == "scatterplot":
        return fn(df, ["cond1"], "pos", **base)
    if name == "boxplot":
        return fn(df, ["cond1"], "muttype", ["Missense", "Silent"], **base)
    if name == "corr_jointplot":
        return fn(df, "cond1", "cond2", **base)
    if name == "corr_heatmap":
        return fn(df, ["cond1", "cond2"], **base)
    if name == "sgrna_enrichment":
        return fn(df, ["cond1"], **base)
    raise ValueError(name)


@pytest.mark.parametrize("name", STATIC_PLOTS)
def test_every_static_plot_returns_its_figure(name, df):
    """Nothing returned a figure before, so no plot could be composed."""
    result = call(name, df)
    assert result is not None, f"{name} returned None"
    fig, _ = result
    assert hasattr(fig, "savefig"), f"{name} did not return a matplotlib figure"


@pytest.mark.parametrize("name", STATIC_PLOTS)
def test_every_static_plot_accepts_a_path(name, df, tmp_path):
    """df_filepath is documented as a path; all of them must honour it."""
    csv = tmp_path / "data.csv"
    df.to_csv(csv, index=False)
    assert call(name, str(csv)) is not None


@pytest.mark.parametrize("name", STATIC_PLOTS)
def test_no_static_plot_mutates_the_caller_frame(name, df):
    """
    Normalizing to negative controls rewrites the comparison columns.

    norm_to_intergenic_ctrls subtracts in place, so a plot that did not copy
    first would silently alter the caller's data.
    """
    before = df.copy()
    call(name, df, neg_ctrl=True, neg_ctrl_col="gene",
         neg_ctrl_conditions=["Control"])
    pd.testing.assert_frame_equal(df, before)


@pytest.mark.parametrize("name", STATIC_PLOTS)
def test_every_static_plot_can_normalize_to_controls(name, df):
    assert call(name, df, neg_ctrl=True, neg_ctrl_col="gene",
                neg_ctrl_conditions=["Control"]) is not None


@pytest.mark.parametrize("name", STATIC_PLOTS)
def test_every_static_plot_can_filter_by_value(name, df):
    assert call(name, df, filter_val=True, val_cols=["score"],
                val_min=2.0) is not None


@pytest.mark.parametrize("name", STATIC_PLOTS)
def test_every_static_plot_can_filter_by_category(name, df):
    assert call(name, df, filter_params=True, params_cols=["muttype"],
                params_conditions=[["Missense", "Silent"]]) is not None


@pytest.mark.parametrize("name", STATIC_PLOTS)
def test_every_static_plot_writes_where_told(name, df, tmp_path):
    out = tmp_path / "nested" / "dir"
    call(name, df, savefig=True, out_dir=str(out), out_name="fig",
         out_type="png")
    written = list(out.glob("*.png"))
    assert written, f"{name} wrote nothing into {out}"


@pytest.mark.parametrize("name", ["scatterplot", "boxplot", "corr_jointplot",
                                  "sgrna_enrichment"])
def test_every_point_plot_can_highlight(name, df):
    assert call(name, df, highlight=True, highlight_col="sgRNA_ID",
                highlight_vals=["sg1", "sg2", "sg3"]) is not None


@pytest.mark.parametrize("name", ["scatterplot", "boxplot", "corr_jointplot",
                                  "sgrna_enrichment"])
def test_every_point_plot_can_colour_by_category(name, df):
    assert call(name, df, include_hue=True, hue_col="muttype") is not None


@pytest.mark.parametrize("name", ["scatterplot", "boxplot"])
@pytest.mark.parametrize("annot_kw", [
    {"annot_top": 5}, {"annot_top": -5}, {"annot_cutoff": 1.0}, {"annot_abs": 5},
])
def test_every_annotating_plot_supports_all_annotation_modes(name, df, annot_kw):
    assert call(name, df, annot=True, annot_label="sgRNA_ID",
                **annot_kw) is not None
