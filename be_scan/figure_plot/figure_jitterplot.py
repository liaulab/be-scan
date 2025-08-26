
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
        ax.set(ylim=axis.ylim, yticks=axis.yticks)

        plt.legend([],[], frameon=False)

        out_path = Path(f"{str(output.path)}_Boxplot.{output.out_type}")
        plt.savefig(out_path, dpi=output.dpi,
                    format=output.out_type, transparent=output.transparent,
                    bbox_inches="tight")
        if output.show: plt.show()
        plt.close()

    if kde:
        # Plot KDEs of Q575R
        fig, ax = plt.subplots(figsize=output.figsize, **output.subplots_kws)

        quartiles =sorted(data['quartile'].unique())[::-1]
        colors = style.colors[::-1]
        for q, c in zip(quartiles, colors):
            print(q, c)
            subset = data[data['quartile'] == q]
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
        ax.set(ylim=axis.ylim, yticks=axis.yticks)

        legend = ax.legend(prop=arial_font6, frameon=False)

        out_path = Path(f"{str(output.path)}_Density.{output.out_type}")
        plt.savefig(out_path, dpi=output.dpi,
                    format=output.out_type, transparent=output.transparent,
                    bbox_inches="tight")
        if output.show: plt.show()
        plt.close()
