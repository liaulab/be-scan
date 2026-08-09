"""
Every figure_plot figure function has an interactive counterpart.

The pair must be interchangeable: same positional arguments, same option
objects, HTML written beside the static file. These tests check the pairing,
the shared styling, and that the interactive versions honour output.save and
output.show the way the static ones do.
"""

import inspect

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

import be_scan.figure_plot as fp
from be_scan.figure_plot.figure_classes import (
    AxisLabelOpts, DomainOpts, OutputOpts, ScatterStyleOpts, resolve_colors,
)
from be_scan.figure_plot import figure_interactive_style as ivs


# static name -> interactive name
PAIRS = {
    "scatterplot_figure": "interactive_scatterplot_figure",
    "correlation_scatterplot_figure": "interactive_correlation_scatterplot_figure",
    "boxplot_figure": "interactive_boxplot_figure",
    "jitterbox_kdeplot_figure": "interactive_jitterbox_kdeplot_figure",
    "lollipop_figure": "interactive_lollipop_figure",
    "lollipop_split_figure": "interactive_lollipop_split_figure",
    "kmeans_clustered_heatmap": "interactive_kmeans_clustered_heatmap",
    "loess_figure": "interactive_loess_figure",
    "pwes_histogram": "interactive_pwes_histogram",
    "pwes_boxplots": "interactive_pwes_boxplots",
    "pwes_scatterplots": "interactive_pwes_scatterplots",
    "pwes_heatmap": "interactive_pwes_heatmap",
    "pwes_clustermap": "interactive_pwes_clustermap",
    "pwes_heatmap_colorbar": "interactive_pwes_heatmap_colorbar",
}

# parameters the interactive versions add and the static ones do not
INTERACTIVE_ONLY = {"scale", "hover_cols"}


@pytest.fixture
def df():
    rng = np.random.default_rng(0)
    n = 200
    return pd.DataFrame({
        "pos": np.arange(1, n + 1),
        "score": rng.normal(size=n),
        "score2": rng.normal(size=n),
        "mut_type": rng.choice(["Missense", "Silent", "Nonsense", "Splice"], n),
        "cat": rng.choice(list("ABC"), n),
        "q": rng.choice([1, 2, 3, 4], n),
        "cluster": rng.integers(1, 5, n),
    })


PAL = {"Missense": "#fb8072", "Silent": "#80b1d3",
       "Nonsense": "#fdb462", "Splice": "#b3de69"}
BPAL = {"A": "#fb8072", "B": "#80b1d3", "C": "#fdb462"}
CPAL = list(PAL.values())


# --------------------------------------------------------------------------- #
# the pairing itself
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("static_name,interactive_name", sorted(PAIRS.items()))
def test_every_figure_function_has_an_interactive_twin(static_name, interactive_name):
    assert hasattr(fp, static_name), f"{static_name} missing"
    assert hasattr(fp, interactive_name), (
        f"{static_name} has no interactive counterpart {interactive_name}"
    )


@pytest.mark.parametrize("static_name,interactive_name", sorted(PAIRS.items()))
def test_the_twin_takes_the_same_arguments(static_name, interactive_name):
    """
    A caller must be able to swap one for the other without editing anything
    else, so every static parameter has to exist on the interactive version
    with the same default.
    """
    static = inspect.signature(getattr(fp, static_name)).parameters
    interactive = inspect.signature(getattr(fp, interactive_name)).parameters

    for name, param in static.items():
        assert name in interactive, (
            f"{interactive_name} is missing '{name}', which {static_name} takes"
        )
        if param.default is not inspect.Parameter.empty:
            assert interactive[name].default == param.default, (
                f"{interactive_name}'s default for '{name}' differs from "
                f"{static_name}'s"
            )

    extra = set(interactive) - set(static) - INTERACTIVE_ONLY
    assert not extra, f"{interactive_name} adds unexpected parameters {sorted(extra)}"


@pytest.mark.parametrize("interactive_name", sorted(PAIRS.values()))
def test_every_twin_offers_scale(interactive_name):
    """scale keeps print-sized figures legible in a browser."""
    params = inspect.signature(getattr(fp, interactive_name)).parameters
    assert "scale" in params
    assert params["scale"].default == 2.0


# --------------------------------------------------------------------------- #
# output behaviour
# --------------------------------------------------------------------------- #


def call_pair(name, df, output, **kw):
    """Invoke an interactive function with the arguments it needs."""
    fn = getattr(fp, name)
    if name == "interactive_scatterplot_figure":
        return fn(df, ["score"], "pos", output=output, **kw)
    if name == "interactive_correlation_scatterplot_figure":
        return fn(df, "score", "score2", list(PAL), hue_col="mut_type",
                  palette=PAL, output=output, **kw)
    if name == "interactive_boxplot_figure":
        return fn(df, "cat", "score", palette=BPAL, output=output, **kw)
    if name == "interactive_jitterbox_kdeplot_figure":
        return fn(df, "q", "score", "cat", BPAL, output=output, **kw)
    if name == "interactive_pwes_histogram":
        return fn(df, CPAL, output=output, **kw)
    if name == "interactive_pwes_boxplots":
        return fn(df, "score", CPAL, output=output, **kw)
    if name == "interactive_pwes_scatterplots":
        return fn(df, "pos", "score", CPAL, output=output, **kw)
    if name == "interactive_kmeans_clustered_heatmap":
        return fn(df[["score", "score2"]], 4, df["cluster"], df["cat"],
                  output=output, **kw)
    raise ValueError(name)


RUNNABLE = [
    "interactive_scatterplot_figure",
    "interactive_correlation_scatterplot_figure",
    "interactive_boxplot_figure",
    "interactive_jitterbox_kdeplot_figure",
    "interactive_pwes_histogram",
    "interactive_pwes_boxplots",
    "interactive_pwes_scatterplots",
    "interactive_kmeans_clustered_heatmap",
]


@pytest.mark.parametrize("name", RUNNABLE)
def test_html_is_written_beside_the_static_figure(name, df, tmp_path):
    """
    The HTML shares the static figure's stem and directory.

    That is the whole point of writing to output.path: a reader finds
    screen.svg and screen.html side by side.
    """
    out = OutputOpts(save=True, show=False, path=str(tmp_path / "nested" / "fig"),
                     out_type="svg")
    call_pair(name, df, out)

    written = list((tmp_path / "nested").glob("*.html"))
    assert written, f"{name} wrote no HTML into {tmp_path/'nested'}"
    for path in written:
        assert path.name.startswith("fig"), (
            f"{name} wrote {path.name}, which does not share the static stem"
        )


@pytest.mark.parametrize("name", RUNNABLE)
def test_save_false_writes_nothing(name, df, tmp_path):
    out = OutputOpts(save=False, show=False, path=str(tmp_path / "fig"))
    call_pair(name, df, out)
    assert not list(tmp_path.glob("*.html")), f"{name} wrote despite save=False"


@pytest.mark.parametrize("name", RUNNABLE)
def test_a_figure_is_returned(name, df, tmp_path):
    out = OutputOpts(save=False, show=False, path=str(tmp_path / "fig"))
    result = call_pair(name, df, out)
    assert result is not None, f"{name} returned None"


def test_html_is_self_contained_by_default(df, tmp_path):
    """
    plotly.js is embedded so a figure opens with no network.

    Large, but a figure emailed to a collaborator has to work on its own.
    """
    out = OutputOpts(save=True, show=False, path=str(tmp_path / "fig"))
    call_pair("interactive_scatterplot_figure", df, out)
    html = (tmp_path / "fig.html").read_text()
    assert "plotly" in html.lower()
    assert len(html) > 1_000_000, "plotly.js does not appear to be embedded"


# --------------------------------------------------------------------------- #
# stylistic consistency with the static figures
# --------------------------------------------------------------------------- #


def test_layout_is_derived_from_the_matplotlib_rcparams():
    """
    Restyling the static figures should restyle the interactive ones.

    The plotly layout reads the live rcParams rather than repeating the
    numbers, so the two cannot drift apart.
    """
    assert ivs.BASE_FONT == matplotlib.rcParams["font.size"]
    assert ivs.BASE_AXES_LW == matplotlib.rcParams["axes.linewidth"]
    assert ivs.BASE_TICK_W == matplotlib.rcParams["xtick.major.width"]
    assert ivs.BASE_TICK_LEN == matplotlib.rcParams["xtick.major.size"]


def test_only_left_and_bottom_spines_are_drawn():
    """The static figures hide the top and right spines; so must these."""
    style = ivs.axis_style(scale=1.0)
    assert style["showline"] is True
    assert style["mirror"] is False, "mirror=True would draw the opposite spine"
    assert style["showgrid"] is False
    assert style["ticks"] == "outside"


def test_figure_size_matches_output_figsize():
    """OutputOpts.figsize is in inches; the pixel size must follow it."""
    out = OutputOpts(figsize=(4.0, 2.2))
    w, h = ivs.figure_size(out, scale=1.0)
    assert w == pytest.approx(4.0 * ivs.PX_PER_INCH)
    assert h == pytest.approx(2.2 * ivs.PX_PER_INCH)

    w2, h2 = ivs.figure_size(out, scale=2.0)
    assert (w2, h2) == pytest.approx((w * 2, h * 2)), "scale must be proportional"


def test_scale_preserves_proportions():
    """
    Doubling scale doubles size and font together.

    If they scaled differently the interactive figure would not look like the
    static one, which is the whole requirement.
    """
    out = OutputOpts(figsize=(4.0, 2.2))
    small = ivs.base_layout(out, scale=1.0)
    large = ivs.base_layout(out, scale=2.0)
    assert large["width"] == small["width"] * 2
    assert large["height"] == small["height"] * 2
    assert large["font"]["size"] == small["font"]["size"] * 2


def test_marker_size_converts_area_to_diameter():
    """
    Matplotlib's `s` is an area in points squared; plotly's size is a
    diameter. Passing one through as the other makes the points far too big.
    """
    assert ivs.marker_size(9, scale=1.0) == pytest.approx(3.0)
    assert ivs.marker_size(25, scale=1.0) == pytest.approx(5.0)
    assert ivs.marker_size(9, scale=2.0) == pytest.approx(6.0)


def test_matplotlib_colour_shorthands_are_translated():
    """
    The style dataclasses use matplotlib conventions, including 'k'.

    Plotly rejects those outright, which is why every colour goes through
    to_css first.
    """
    assert ivs.to_css("k") == "#000000"
    assert ivs.to_css("r") == "#ff0000"
    assert ivs.to_css("#fb8072") == "#fb8072"
    assert ivs.to_css("lightblue") == "#add8e6"
    assert ivs.to_css(None) == "black"
    assert ivs.to_css("not-a-colour") == "black"


def test_line_style_translates_to_plotly_dashes():
    line = ivs.line_from_kws({"color": "k", "ls": "--", "lw": 1}, scale=2.0)
    assert line["color"] == "#000000"
    assert line["dash"] == "dash"
    assert line["width"] == 2.0
    assert ivs.line_from_kws({"ls": ":"})["dash"] == "dot"
    assert ivs.line_from_kws({"ls": "-"})["dash"] == "solid"


def test_categories_get_the_same_colour_in_both_renderings():
    """
    A reader switching from the static figure to the interactive one must see
    a category keep its colour. Both call the same resolver on sorted values.
    """
    values = ["Silent", "Missense", "Nonsense"]
    style = ScatterStyleOpts()
    first = resolve_colors(values, style)
    second = resolve_colors(list(reversed(values)), style)
    assert first == second, "colour assignment depends on input order"

    explicit = resolve_colors(values, style, PAL)
    assert explicit["Missense"] == PAL["Missense"]


def test_palette_string_is_not_iterated_character_by_character():
    """
    ScatterStyleOpts.palette defaults to the string "colorblind".

    Zipping values against it yields the characters 'c', 'o', 'l', which is
    what the static scatterplot used to do.
    """
    colors = resolve_colors(["a", "b", "c"], ScatterStyleOpts())
    assert all(c.startswith("#") for c in colors.values()), colors
    assert len(set(colors.values())) == 3


@pytest.mark.parametrize("name", RUNNABLE)
def test_no_twin_mutates_the_caller_frame(name, df, tmp_path):
    before = df.copy()
    out = OutputOpts(save=False, show=False, path=str(tmp_path / "fig"))
    call_pair(name, df, out)
    pd.testing.assert_frame_equal(df, before)
