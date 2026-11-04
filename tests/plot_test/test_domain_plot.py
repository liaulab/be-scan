"""
domain_plot_figure draws a residue-proportional domain map.

The properties that matter: one residue is the same width everywhere, colours
fall back to the default unless they are real hex codes, and no two labels
overlap whatever the angle.
"""

import itertools
from dataclasses import replace

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import pytest
from matplotlib.patches import FancyBboxPatch

import be_scan.figure_plot as fp
from be_scan.figure_plot.figure_domain_plot import is_hex_color, spread_labels


SMARCA4 = [
    ("IDR1", 1, 175), ("QLQ", 176, 207), ("IDR1", 208, 338),
    ("pre-HSA", 339, 422), ("HSA", 423, 508), ("Unk.", 509, 566),
    ("(unannotated)", 567, 593), ("BRK", 594, 630), ("(unannotated)", 631, 653),
    ("pre-Helicase", 654, 710), ("Helicase-N", 711, 955),
    ("Helicase-C", 956, 1229), ("SnAC", 1230, 1355), ("IDR2", 1356, 1381),
    ("Bromodomain", 1382, 1510), ("IDR2", 1511, 1590),
]
COLORS = ["#fb8072", "#80b1d3"]


@pytest.fixture
def domains():
    return [dict(description=d, start=s, end=e, color=COLORS[i % 2])
            for i, (d, s, e) in enumerate(SMARCA4)]


@pytest.fixture
def out(tmp_path):
    return replace(fp.DOMAINPLOTOUTPUT, save=False, show=False,
                   path=str(tmp_path / "dom"))


def boxes(ax):
    return [p for p in ax.patches if isinstance(p, FancyBboxPatch)]


def test_default_size_is_scatterplot_width_and_a_tenth_tall():
    w, h = fp.DOMAINPLOTOUTPUT.figsize
    assert w == fp.OUTPUT.figsize[0]
    assert h == pytest.approx(w / 10)


def test_one_box_per_domain_proportional_to_its_length(domains, out):
    fig, ax = fp.domain_plot_figure(domains, 1590, output=out)
    drawn = boxes(ax)
    assert len(drawn) == len(domains)
    for patch, d in zip(drawn, domains):
        assert patch.get_x() == d["start"] - 1
        assert patch.get_width() == d["end"] - d["start"] + 1


def test_axes_span_the_whole_figure_width(domains, out):
    fig, ax = fp.domain_plot_figure(domains, 1590, output=out)
    assert ax.get_position().x0 == 0 and ax.get_position().x1 == 1
    assert ax.get_xlim() == (0, 1590)


@pytest.mark.parametrize("color,expected", [
    ("#fb8072", "#fb8072"), ("#ABC", "#aabbcc"),
    ("red", "#b7b7b7"), ("fb8072", "#b7b7b7"), (None, "#b7b7b7"),
    ((1, 0, 0), "#b7b7b7"),
])
def test_only_hex_codes_override_the_default_colour(color, expected, out):
    fig, ax = fp.domain_plot_figure(
        [dict(description="A", start=1, end=10, color=color)], 100, output=out)
    assert mcolors.to_hex(boxes(ax)[0].get_facecolor()) == expected


def test_defaults_black_half_point_outline_and_baseline(out):
    fig, ax = fp.domain_plot_figure([dict(start=1, end=10)], 100, output=out)
    box = boxes(ax)[0]
    assert mcolors.to_hex(box.get_edgecolor()) == "#000000"
    assert box.get_linewidth() == 0.5
    base = ax.lines[0]
    assert mcolors.to_hex(base.get_color()) == "#000000"
    assert base.get_linewidth() == 0.5
    assert list(base.get_xdata()) == [0, 100]


@pytest.mark.parametrize("angle", [0, 30, 45, 90, -45])
def test_labels_never_overlap(domains, out, angle):
    fig, ax = fp.domain_plot_figure(domains, 1590, label_angle=angle, output=out)
    renderer = fig.canvas.get_renderer()
    texts = ax.texts
    assert [t.get_text() for t in texts] == [d for d, _, _ in SMARCA4]
    assert all(t.get_rotation() == pytest.approx(angle % 360) for t in texts)
    if angle % 90 == 0:
        # axis-aligned: the rendered boxes themselves must be disjoint
        for a, b in itertools.combinations(texts, 2):
            assert not a.get_window_extent(renderer).overlaps(
                b.get_window_extent(renderer))
    # labels stay in domain order left to right
    xs = [t.get_position()[0] for t in texts]
    assert xs == sorted(xs)


def test_labels_use_the_figure_font(domains, out):
    fig, ax = fp.domain_plot_figure(domains, 1590, output=out)
    assert all(t.get_fontsize() == 6 for t in ax.texts)


def test_saved_svg_keeps_labels_as_text(domains, tmp_path):
    out = replace(fp.DOMAINPLOTOUTPUT, save=True, show=False,
                  path=str(tmp_path / "dom"))
    fp.domain_plot_figure(domains, 1590, label_angle=45, output=out)
    svg = (tmp_path / "dom.svg").read_text()
    assert "Bromodomain" in svg and "Helicase-N" in svg


@pytest.mark.parametrize("bad", [
    [dict(start=10, end=5)], [dict(start=1)],
])
def test_bad_domains_are_rejected(bad, out):
    with pytest.raises(ValueError):
        fp.domain_plot_figure(bad, 100, output=out)


def test_bad_gene_length_is_rejected(out):
    with pytest.raises(ValueError):
        fp.domain_plot_figure([dict(start=1, end=5)], 0, output=out)


def test_is_hex_color():
    assert is_hex_color("#fff") and is_hex_color("#B7B7B7") and is_hex_color("#b7b7b7ff")
    assert not is_hex_color("gray") and not is_hex_color("#ggg") and not is_hex_color(None)


def test_spread_labels_leaves_separated_labels_alone():
    assert spread_labels([0, 100], [10, 10]) == [0, 100]


def test_spread_labels_centres_a_colliding_group_on_its_anchors():
    centres = spread_labels([50, 50, 50], [10, 10, 10], gap=2)
    assert centres == pytest.approx([38, 50, 62])
