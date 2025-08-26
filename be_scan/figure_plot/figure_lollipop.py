
from be_scan.figure_plot.figure_classes import *

def lollipop_figure(
    # POSITIVE #
    df_pos,
    pos_xcol, pos_ycol, pos_threshold,
    pos_category_col, pos_colors,

    # NEGATIVE #
    df_neg,
    neg_xcol, neg_ycol, neg_threshold,
    neg_category_col, neg_colors,

    axis: AxisLabelOpts = AXIS,
    domain: DomainOpts = DOMAIN,
    style: LollipopStyleOpts = LOLLIPOPSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    fig, ax = plt.subplots(figsize=output.figsize)

    for spine in ax.spines.values():
        spine.set_linewidth(axis.linewidth)

    # PLOTTING HITS ONLY #
    pos_hits = df_pos.loc[df_pos[pos_ycol] > pos_threshold].copy()
    neg_hits = df_neg.loc[df_neg[neg_ycol] < neg_threshold].copy()
    print('Positive Hit Count:', len(pos_hits))
    print('Negative Hit Count:', len(neg_hits))

    # DRAW DOMAINS #
    for d in domain.dom_setting:
        ax.axvspan(d['start'], d['end'], facecolor=d['color'],
                   alpha=domain.default_alpha, )

    # POSITIVE #
    pos_unique = pos_hits[pos_category_col].unique()
    for category, color in zip(pos_unique, pos_colors):

        subset = pos_hits[pos_hits[pos_category_col] == category]
        mark, stem, base = ax.stem(subset[pos_xcol], subset[pos_ycol],
                                   label = f"Pos {category}", basefmt=" ")

        plt.setp(mark, markerfacecolor = color, **style.marker_kws)
        plt.setp(stem, **style.stem_kws)
        plt.setp(base, **style.base_kws)
        mark.set_rasterized(style.marker_rasterized)
        stem.set_rasterized(style.stem_rasterized)

    # NEGATIVE #
    neg_unique = neg_hits[neg_category_col].unique()
    for category, color in zip(neg_unique, neg_colors):

        subset = neg_hits[neg_hits[neg_category_col] == category]
        mark, stem, base = ax.stem(subset[neg_xcol], subset[neg_ycol],
                                   label = f"Neg {category}", basefmt=" ")

        plt.setp(mark, markerfacecolor = color, **style.marker_kws)
        plt.setp(stem, **style.stem_kws)
        plt.setp(base, **style.base_kws)
        mark.set_rasterized(style.marker_rasterized)
        stem.set_rasterized(style.stem_rasterized)

    # SPINES #
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    ax.set_xticks(axis.xticks)
    ax.set_xlim(axis.xlim)
    ax.set_ylim(axis.ylim)
    ax.set_yticks(axis.yticks)
    for label in ax.get_xticklabels() + ax.get_yticklabels(): label.set_fontproperties(arial_font6)

    # AXIS LABELS AND LEGEND #
    ax.axhline(y=0, **style.line_kws)
    ax.set_title(axis.title, fontproperties=arial_font6)
    ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    legend = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    handles, labels = ax.get_legend_handles_labels()
    legend.remove()

    # SAVE #
    if output.save:
        out_path = Path(f"{str(output.path)}.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent
                    )

    # SHOW #
    if output.show: plt.show()
    plt.close(fig)

    return fig, legend, handles, labels

def lollipop_split_figure(
    # POSITIVE #
    df_pos,
    pos_xcol, pos_ycol, pos_threshold,
    pos_category_col, pos_colors,
    ylim_top, yticks_top,

    # NEGATIVE #
    df_neg,
    neg_xcol, neg_ycol, neg_threshold,
    neg_category_col, neg_colors,
    ylim_bot, yticks_bot,

    height_ratios = [1, 0.001, 1],

    axis: AxisLabelOpts = AXIS,
    domain: DomainOpts = DOMAIN,
    style: LollipopStyleOpts = LOLLIPOPSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    fig = plt.figure(figsize=output.figsize)
    gs = gridspec.GridSpec(3, 1, height_ratios=height_ratios)
    ax_top = fig.add_subplot(gs[0])
    ax_bot = fig.add_subplot(gs[2])

    # PLOTTING HITS ONLY #
    pos_hits = df_pos.loc[df_pos[pos_ycol] > pos_threshold].copy()
    neg_hits = df_neg.loc[df_neg[neg_ycol] < neg_threshold].copy()
    print('Positive Hit Count:', len(pos_hits))
    print('Negative Hit Count:', len(neg_hits))

    # DRAW DOMAINS #
    for d in domain.dom_setting:
        for ax in (ax_top, ax_bot):
            ax.axvspan(d['start'], d['end'], facecolor=d['color'],
                       alpha=domain.default_alpha)

    # POSITIVE #
    pos_unique = pos_hits[pos_category_col].unique()
    for category, color in zip(pos_unique, pos_colors):

        subset = pos_hits[pos_hits[pos_category_col] == category]
        mark, stem, base = ax_top.stem(subset[pos_xcol], subset[pos_ycol],
                                   label = f"Pos {category}", basefmt=" ")

        plt.setp(mark, markerfacecolor = color, **style.marker_kws)
        plt.setp(stem, **style.stem_kws)
        plt.setp(base, **style.base_kws)
        mark.set_rasterized(style.marker_rasterized)
        stem.set_rasterized(style.stem_rasterized)

    # NEGATIVE #
    neg_unique = neg_hits[neg_category_col].unique()
    for category, color in zip(neg_unique, neg_colors):

        subset = neg_hits[neg_hits[neg_category_col] == category]
        mark, stem, base = ax_bot.stem(subset[neg_xcol], subset[neg_ycol],
                                   label = f"Neg {category}", basefmt=" ")

        plt.setp(mark, markerfacecolor = color, **style.marker_kws)
        plt.setp(stem, **style.stem_kws)
        plt.setp(base, **style.base_kws)
        mark.set_rasterized(style.marker_rasterized)
        stem.set_rasterized(style.stem_rasterized)

    # SPINES #
    for ax in (ax_top, ax_bot):
        for spine in ax.spines.values():
            spine.set_linewidth(axis.linewidth)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_position(('outward', 3))
        ax.set_xticks(axis.xticks)
        ax.set_xlim(axis.xlim)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontproperties(arial_font6)

    # AXIS LABELS AND LEGEND #
    ax_top.spines['bottom'].set_visible(False)
    ax_top.set_xticklabels([])
    ax_top.axhline(y=0, **style.line_kws)
    ax_bot.axhline(y=0, **style.line_kws)
    ax_top.set_ylim(ylim_top[0], ylim_top[1]) ###
    ax_bot.set_ylim(ylim_bot[0], ylim_bot[1]) ###
    ax_top.set_yticks(yticks_top) ###
    ax_bot.set_yticks(yticks_bot) ###
    ax_bot.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    ax_top.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    ax_top.set_title(axis.title, fontproperties=arial_font6)

    # LEGEND #
    handles_top, labels_top = ax_top.get_legend_handles_labels()
    handles_bot, labels_bot = ax_bot.get_legend_handles_labels()
    handles = handles_top + handles_bot
    labels = labels_top + labels_bot
    legend = fig.legend(handles, labels, bbox_to_anchor=(1.05, 1), loc='upper left')
    legend.remove()

    # SAVE #
    if output.save:
        out_path = Path(f"{str(output.path)}.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent
                    )

    # SHOW #
    if output.show: plt.show()
    plt.close(fig)

    return fig, legend, handles, labels

