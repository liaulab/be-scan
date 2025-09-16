
from be_scan.figure_plot.figure_classes import *
from be_scan.figure_plot.figure_loess_helpers import *

def loess_figure(
    results: dict,
    gene : str = '',
    max_cap : int = 4,
    pval_thresh : float = 0.05,
    *,
    axis_pval: AxisLabelOpts = AXIS,
    axis_regr: AxisLabelOpts = AXIS,
    style: LoessStyleOpts = LOESSSTYLE,
    domain: DomainOpts = DOMAIN,
    output_pval: OutputOpts = OUTPUT,
    output_regr: OutputOpts = OUTPUT,
    ):

    for comp, res in results.items():
        df_pvals = res['pvals']
        df_loess = res['loess']

        # --- P-VALUE PLOT --- #

        fig_pval, ax_pval = plt.subplots(figsize=output_pval.figsize, constrained_layout=False)
        df_plotting = df_pvals[['corr_pval']].copy()
        df_plotting['-log10'] = -1 * np.log10(df_plotting['corr_pval'] + 10**-max_cap)
        threshold = -1 * np.log10(pval_thresh + 10**-max_cap) # -log(p=0.05) #

        # X Y VALUES #
        xvals, yvals = df_plotting.index.to_numpy(), df_plotting['-log10'].to_numpy()
        xvals, yvals = increase_resolution(xvals, yvals, threshold)
        x_listoflists, y_listoflists = split_on_value(xvals, yvals, threshold)

        # COLOR DIFFERENTLY ABOVE AND BELOW THRESHOLD #
        for x_list, y_list in zip(x_listoflists, y_listoflists):
            x_array, y_array = np.array(x_list, dtype=float), np.array(y_list, dtype=float)

            # PLOT BELOW THRESHOLD #
            if max(y_array) <= threshold: color = style.below_color
            # PLOT ABOVE THRESHOLD #
            if max(y_array) > threshold: color = style.above_color

            ax_pval.plot(x_array, y_array, color=color, **style.pval_plot_kws)
        ax_pval.axhline(y=threshold, **style.line_kws)
        loess_figure_style(fig_pval, ax_pval, axis_pval, domain, output_pval)

        # --- LOESS REGRESSION PLOT --- #

        fig_regr, ax_regr = plt.subplots(figsize=output_regr.figsize, constrained_layout=False)
        df_plotting = df_loess[['y_loess']].copy()
        # PLOT #
        ax_regr.plot(df_plotting.index, df_plotting['y_loess'], **style.regr_plot_kws)
        ax_regr.axhline(y=0, **style.line_kws)
        loess_figure_style(fig_regr, ax_regr, axis_regr, domain, output_regr)

    return fig_pval, ax_pval, fig_regr, ax_regr

def loess_figure_style(
    fig, ax, axis, domain, output
    ):

    # LABEL SETTINGS #
    ax.set_title(axis.title, fontproperties=arial_font6)
    ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)

    # AXES SETTINGS #
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    ax.set_xlim(axis.xlim)
    ax.set_ylim(axis.ylim)
    if axis.xticks: ax.set_xticks(axis.xticks)
    if axis.yticks: ax.set_yticks(axis.yticks)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(arial_font6)
    for col in ax.collections:
        if hasattr(col, 'set_linewidth'): col.set_linewidth(axis.linewidth)
    for spine in ax.spines.values(): spine.set_linewidth(axis.linewidth)
    for line in ax.get_lines(): line.set_linewidth(axis.linewidth)

    # PLOT DOMAINS #
    for dom in domain.dom_setting:
        ax.axvspan(
            dom["start"], dom["end"],
            facecolor=dom.get("color", domain.default_color),
            edgecolor="none",
            alpha=dom.get("alpha", domain.default_alpha) )

    # SAVE AND CLOSE #
    fig.tight_layout()
    if output.save:
        # scale_figure_elements(fig1, arial_font6)
        out_path = Path(f"{str(output.path)}.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                          dpi=output.dpi, transparent=output.transparent
                          )
    if output: plt.show(fig)
