"""
Author: Calvin XiaoYang Hu

{Description: a residue-proportional linear map of a protein's domains}

Every residue gets the same horizontal width, so a domain's on-page length is
proportional to its length in residues. The axes fill the whole figure width,
which makes the scale exactly figsize[0] / gene_len inches per residue and lets
the bar line up with a scatterplot or lollipop of the same width drawn above
or below it.

    gene_domains = [
        {'description': 'IDR1', 'start': 1,   'end': 175, 'color': '#fb8072'},
        {'description': 'QLQ',  'start': 176, 'end': 207, 'color': '#80b1d3'},
        ...
    ]
    domain_plot_figure(gene_domains, gene_len=1590, label_angle=45)
"""

import re

import matplotlib.transforms as mtransforms
from matplotlib.patches import FancyBboxPatch

from be_scan.figure_plot.figure_classes import *

_HEX_COLOR = re.compile(
    r'^#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{4}|[0-9a-fA-F]{6}|[0-9a-fA-F]{8})$')


def is_hex_color(color):
    """True only for a '#rgb', '#rgba', '#rrggbb' or '#rrggbbaa' string."""
    return isinstance(color, str) and bool(_HEX_COLOR.match(color.strip()))


def domain_fill(dom, default_color):
    """The fill for one domain: its own colour if it is a valid hex code."""
    color = dom.get('color')
    return color.strip() if is_hex_color(color) else default_color


def spread_labels(anchors, widths, gap=0.0):
    """
    Place 1D labels as close to their anchors as possible without overlap.

    anchors are the preferred centres, widths the room each label needs (both
    in the same units). Overlapping labels are merged into a block centred on
    the mean of its members' anchors, and blocks keep merging until none
    overlap. Returns the centres in the original order.
    """
    n = len(anchors)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: anchors[i])

    # block: (members, offset of each member's centre from the block's left
    #         edge, total width, left edge)
    blocks = []
    for i in order:
        w = widths[i]
        members, offsets, total, left = [i], [w / 2], w, anchors[i] - w / 2
        while blocks:
            p_members, p_offsets, p_total, p_left = blocks[-1]
            if p_left + p_total + gap <= left:
                break
            blocks.pop()
            shift = p_total + gap
            members = p_members + members
            offsets = p_offsets + [o + shift for o in offsets]
            total = p_total + gap + total
            left = float(np.mean([anchors[m] - o for m, o in zip(members, offsets)]))
        blocks.append((members, offsets, total, left))

    centres = [0.0] * n
    for members, offsets, _, left in blocks:
        for m, o in zip(members, offsets):
            centres[m] = left + o
    return centres


@dataclass(frozen = True)
class DomainPlotStyleOpts:
    # used for any domain whose 'color' is missing or not a valid hex code
    default_color: str = "#b7b7b7"
    rect_kws: Dict[str, Any] = field(default_factory=lambda:
        {"edgecolor": "black", "linewidth": 0.5})
    baseline_kws: Dict[str, Any] = field(default_factory=lambda:
        {"color": "black", "linewidth": 0.5})
    # fraction of the figure height the domain boxes occupy, centred on the baseline
    bar_height: float = 0.6
    # corner radius as a fraction of the box height; 0 for square corners
    corner_radius: float = 0.18
    # points between the bottom of the boxes and the top of the labels
    label_pad: float = 1.5
    # points of empty space kept between neighbouring labels
    label_gap: float = 1.5


DOMAINPLOTSTYLE = DomainPlotStyleOpts()

# the scatterplot width, with a height a tenth of it
DOMAINPLOTOUTPUT = replace(OUTPUT, figsize=(OUTPUT.figsize[0], OUTPUT.figsize[0] / 10),
                           path="domain_plot")


def domain_plot_figure(
    gene_domains: List[Dict[str, Any]],
    gene_len: int,
    *,
    label_angle: float = 0,
    style: DomainPlotStyleOpts = DOMAINPLOTSTYLE,
    output: OutputOpts = DOMAINPLOTOUTPUT,
):
    """
    Draw a residue-proportional domain map with labels underneath.

    Parameters
    ----------
    gene_domains : list of dict
        One dict per domain with 'start' and 'end' (1-based, inclusive), and
        optionally 'description' (the label) and 'color'. A 'color' that is
        not a valid hex code falls back to style.default_color.
    gene_len : int
        Length of the whole protein, in residues. Sets the baseline and the
        horizontal scale.
    label_angle : float
        Rotation of the labels in degrees, counter-clockwise. 0 is horizontal;
        45 or 90 fit more labels under a narrow figure.
    style : DomainPlotStyleOpts
    output : OutputOpts
        Defaults to the scatterplot width and a tenth of it in height. Saved
        with a tight bounding box so the labels below the bar are kept.

    Returns
    -------
    fig, ax
    """
    if gene_len is None or gene_len <= 0:
        raise ValueError("gene_len must be a positive number of residues")
    for dom in gene_domains:
        if 'start' not in dom or 'end' not in dom:
            raise ValueError(f"domain is missing 'start' or 'end': {dom}")
        if dom['end'] < dom['start']:
            raise ValueError(f"domain ends before it starts: {dom}")

    fig = plt.figure(figsize=output.figsize, **output.subplots_kws)
    # axes fill the figure, so one residue is exactly width / gene_len inches
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, gene_len)
    ax.set_ylim(0, 1)
    ax.axis('off')

    # BASELINE #
    ax.plot([0, gene_len], [0.5, 0.5], solid_capstyle='butt', zorder=1,
            clip_on=False, **style.baseline_kws)

    # DOMAINS #
    # residue r spans [r-1, r], so residues 1..gene_len fill 0..gene_len
    bottom = 0.5 - style.bar_height / 2
    width_pt, height_pt = (s * 72 for s in output.figsize)
    x_per_pt, y_per_pt = gene_len / width_pt, 1 / height_pt
    radius_pt = style.corner_radius * style.bar_height * height_pt
    for dom in gene_domains:
        start = max(dom['start'] - 1, 0)
        end = min(dom['end'], gene_len)
        # rounding is given in x units; mutation_aspect keeps corners circular
        # although x and y are on very different data scales
        r_pt = min(radius_pt, (end - start) / x_per_pt / 2,
                   style.bar_height * height_pt / 2)
        ax.add_patch(FancyBboxPatch(
            (start, bottom), end - start, style.bar_height,
            boxstyle=f"round,pad=0,rounding_size={r_pt * x_per_pt}",
            mutation_aspect=y_per_pt / x_per_pt,
            facecolor=domain_fill(dom, style.default_color),
            zorder=2, clip_on=False, **style.rect_kws))

    # LABELS #
    labelled = [d for d in gene_domains if d.get('description')]
    if labelled:
        _draw_labels(fig, ax, labelled, gene_len, bottom, label_angle, style)

    # tight bbox: the labels hang below the figure's nominal height
    finish_figure(fig, output, bbox_inches='tight', pad_inches=0.02)
    return fig, ax


def _draw_labels(fig, ax, domains, gene_len, bar_bottom, angle, style):
    """Labels below the boxes, nudged sideways so none overlap."""
    angle = float(angle)
    rad = np.deg2rad(angle)
    sin, cos = abs(np.sin(rad)), abs(np.cos(rad))
    rotated = sin > 1e-6

    # Horizontal labels sit centred under their domain. Rotated ones hang from
    # their end so the text reads down and away from the bar.
    if rotated:
        ha, va = ('right' if np.sin(rad) > 0 else 'left'), 'center'
    else:
        ha, va = 'center', 'top'

    anchors = [(max(d['start'] - 1, 0) + min(d['end'], gene_len)) / 2 for d in domains]
    below = mtransforms.offset_copy(ax.transData, fig=fig, x=0, y=-style.label_pad,
                                    units='points')
    texts = [ax.text(x, bar_bottom, d['description'], transform=below,
                     fontproperties=arial_font6, rotation=0,
                     rotation_mode='anchor', ha=ha, va=va, clip_on=False)
             for x, d in zip(anchors, domains)]

    # measure unrotated text to know how much horizontal room each one needs
    renderer = fig.canvas.get_renderer()
    px_per_unit = ax.get_window_extent(renderer).width / gene_len
    widths = []
    for t in texts:
        bbox = t.get_window_extent(renderer)
        if not rotated:
            need = bbox.width
        else:
            # parallel strips hanging from the same height clear each other
            # once they are one text height apart perpendicular to the text
            need = bbox.height / sin
            if cos > 1e-6:
                need = min(need, bbox.width / cos)
        widths.append(need / px_per_unit)
        t.set_rotation(angle)

    gap = style.label_gap * fig.dpi / 72 / px_per_unit
    for t, x in zip(texts, spread_labels(anchors, widths, gap=gap)):
        t.set_x(x)
