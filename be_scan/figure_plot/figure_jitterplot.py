
from be_scan.figure_plot.figure_classes import *

def jitterbox_kdeplot_figure(
    data,
    xcol, ycol,
    huecol, color_dict,
    jitterbox = True, kde = True, threshold = None,

    axis: AxisLabelOpts = AXIS,
    style: JitterStyleOpts = JITTERPLOTSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    data = load_frame(data)
    figures = {}

    if jitterbox:
        fig, ax = plt.subplots(figsize=output.figsize, **output.subplots_kws)

        # Plot by WT Dropout
        sns.boxplot(
            data=data,
            x=xcol, y=ycol,
            **style.box_kws,
            boxprops=style.boxprops,
            whiskerprops=style.whiskerprops,
            capprops=style.capprops,
            medianprops=style.medianprops,
        )

        if threshold:
            sns.stripplot(
                data=data[data[ycol] <= threshold],
                x=xcol, y=ycol,
                hue=huecol, palette=color_dict,
                **style.stripplot_kws,
            )
            sns.stripplot(
                data=data[data[ycol] > threshold],
                x=xcol, y=ycol,
                hue=huecol, palette=color_dict,
                **style.stripplot_hits_kws,
            )
        else:
            sns.stripplot(
                data=data,
                x=xcol, y=ycol,
                hue=huecol, palette=color_dict,
                **style.stripplot_kws,
            )

        for spine in ax.spines.values():
            spine.set_linewidth(style.linewidth)

        ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
        ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
        ax.set_title(axis.title, fontproperties=arial_font6)
        plt.xticks(fontproperties=arial_font6)
        plt.yticks(fontproperties=arial_font6)
        apply_axis_limits(ax, axis)

        plt.legend([],[], frameon=False)

        # this saved unconditionally, so output.save had no effect
        finish_figure(fig, output, suffix="_Boxplot", bbox_inches="tight")
        figures["jitterbox"] = (fig, ax)

    if kde:
        # Plot KDEs of Q575R
        fig, ax = plt.subplots(figsize=output.figsize, **output.subplots_kws)

        quartiles =sorted(data[xcol].unique())[::-1]
        colors = style.colors[::-1]
        for q, c in zip(quartiles, colors):
            subset = data[data[xcol] == q]
            sns.kdeplot(
                y=subset[ycol], label=f'Q{q}',
                color=c,
                **style.kdeplot_kws,
                )

        for spine in ax.spines.values():
            spine.set_linewidth(style.linewidth)

        ax.set_xlabel('Density', fontproperties=arial_font6)
        ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
        ax.set_title(axis.title, fontproperties=arial_font6)
        plt.xticks(fontproperties=arial_font6)
        plt.yticks(fontproperties=arial_font6)
        apply_axis_limits(ax, axis)

        legend = ax.legend(prop=arial_font6, frameon=False)

        finish_figure(fig, output, suffix="_Density", bbox_inches="tight")
        figures["kde"] = (fig, ax)

    return figures
