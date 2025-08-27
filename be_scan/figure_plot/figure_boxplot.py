
from be_scan.figure_plot.figure_classes import *

def boxplot_figure(
    df,
    x_col, y_col,
    palette,
    lab_order=None,

    axis: AxisLabelOpts = AXIS,
    style: BoxStyleOpts = BOXPLOTSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    fig, ax = plt.subplots(figsize=output.figsize,
                           **output.subplots_kws)

    # HORIZONTAL LINE #
    plt.axhline(0, **style.axhline_kws)

    # BOXPLOT LINE #
    sns.boxplot(
        data = df,
        ax = ax,
        x = x_col, y = y_col,
        order = lab_order,
        palette=palette, ###

        **style.box_kws,
        boxprops=style.boxprops,
        whiskerprops=style.whiskerprops,
        capprops=style.capprops,
        medianprops=style.medianprops,
    )

    ax.set(ylim=axis.ylim, yticks=axis.yticks)
    ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    ax.set_title(axis.title, fontproperties=arial_font6)

    # Tick labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=style.xlabel_rotation, ha=style.xlabel_ha)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(arial_font6)

    out_path = Path(f"{str(output.path)}.{output.out_type}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=output.dpi,
                format=output.out_type, transparent=output.transparent,
                bbox_inches="tight")
    if output.show: plt.show()

    return fig, ax
