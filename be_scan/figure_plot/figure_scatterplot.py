
from be_scan.figure_plot.figure_classes import *

def scatterplot_figure(
    df_data: Union[pd.DataFrame, str],
    y_col_list: List[str],
    x_col: str,
    *,
    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    domain: DomainOpts = DOMAIN,
    legend: LegendOpts = LEGEND,
    output: OutputOpts = OUTPUT,
    # neg_ctrl: Optional[NegCtrlOpts] = None,
):

    # LOAD DATA OR COPY DF #
    if isinstance(df_data, (str, Path)):
        df_data = pd.read_csv(df_data)
    else: df_data = df_data.copy()

    # MAKE SURE ALL NEEDED COLUMNS ARE PRESENT #
    needed = set(y_col_list + [x_col])
    for col in (style.hue_col, style.marker_col, domain.gene_col):
        if col: needed.add(col)
    missing = needed.difference(df_data.columns)
    if missing: raise ValueError(f"Missing column(s): {sorted(missing)}")

    # FILTER X-VALUES, MIN AND MAX INCLUSIVE #
    if axis.xwindow: df_data = df_data[df_data[x_col].between(*axis.xwindow)]

    # # Negative‑control normalisation
    # ctrl_stats = None
    # if neg_ctrl:
    #     assert(neg_ctrl.neg_ctrl_col in df_data.columns.tolist())
    #     _, ctrl_stats, avg_dict = calc_neg_ctrls(df_data, y_col_list, neg_ctrl.neg_ctrl_col, neg_ctrl.neg_ctrl_conditions)
    #     df_data = norm_to_intergenic_ctrls(df_data, y_col_list, avg_dict)

    # SETUP FIGURE AND AXES #
    n = len(y_col_list)
    fig, axes = plt.subplots(nrows=n, ncols=1,
                             figsize=(output.figsize[0], output.figsize[1]*n),
                             **output.subplots_kws)
    axes = axes if n > 1 else [axes]

    # INCORPORATE COLOR MAP IF ONE IS PROVIDED #
    pal_dict = {}
    if style.hue_col and not style.color_map:
        hue_vals = sorted(df_data[style.hue_col].unique())
        pal_dict = {key: value for key, value in zip(hue_vals, style.palette)}

    # ----- SCATTERPLOT ----- #

    for ax, ycol in zip(axes, y_col_list):

        # PLOT DOMAINS #
        for dom in domain.dom_setting:
            ax.axvspan(
                dom["start"], dom["end"],
                facecolor=dom.get("color", domain.default_color),
                edgecolor="none",
                alpha=dom.get("alpha", domain.default_alpha) )

        # CREATE A MASK FOR TRANSPARENCY BELOW AN ABSOLUTE THRESHOLD #
        transp_mask = df_data[ycol].abs() >= style.transparency_threshold
        layers = [
            (df_data[~transp_mask], style.transparent_kws),
            (df_data[transp_mask],  style.opaque_kws)
        ]

        # PLOT TRANSPARENT AND OPAQUE LAYERS #
        for layer_df, layer_style in layers:
            if layer_df.empty: continue
            base_kws = {**layer_style}
            plot_kws = {
                "x": x_col, "y": ycol, "data": layer_df
                }
            # DEFINE MARKER SETTINGS #
            plot_kws["style"] = style.marker_col
            plot_kws["markers"] = list(style.marker_map.values())
            # DEFINE STYLE SETTINGS #
            if style.hue_col:
                plot_kws["hue"] = style.hue_col
                if style.color_map: plot_kws["palette"] = style.color_map
                else: plot_kws["palette"] = pal_dict
            # PLOT #
            sns.scatterplot(ax=ax, **plot_kws, **base_kws)

        # # If negative control is given, draw in 2 SD guide lines
        # if ctrl_stats and neg_ctrl:
        #     sd = next((s for name, _, s in ctrl_stats if name == ycol), None)
        #     if sd is not None:
        #         ax.axhline(2*sd, **style.ctrl_line_kws)
        #         ax.axhline(-2*sd, **style.ctrl_line_kws)

        if axis.title: ax.set_title(axis.title, fontproperties=arial_font6)

        # SET LIMITS #
        if axis.xlim: ax.set_xlim(*axis.xlim)
        if axis.ylim: ax.set_ylim(*axis.ylim)

        # SET X-TICKS #
        if axis.xticks: ax.set_xticks(axis.xticks)
        elif axis.protein_len: ax.set_xticks([1, axis.protein_len])
        # SET Y-TICKS #
        if axis.yticks: ax.set_yticks(axis.yticks)
        # TICK LABELS #
        for label in ax.get_xticklabels() + ax.get_yticklabels(): label.set_fontproperties(arial_font6)

        # AXES LABELS #
        ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
        ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)

        # LINE SETTINGS #
        for spine in ax.spines.values(): spine.set_linewidth(axis.linewidth)
        for line in ax.get_lines(): line.set_linewidth(axis.linewidth)

        # SPINES AND LEGEND SETTINGS #
        ax.spines['left'].set_position(('outward', 3))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig_leg = ax.get_legend()
        if fig_leg is not None: fig_leg.remove()

    # ----- LEGEND ----- #

    handles, labels = [], []

    # MARKERS #
    if style.marker_col and style.marker_map:
        for val, mark in style.marker_map.items():
            handles.append(plt.Line2D([0], [0], marker=mark, linestyle="", color="black", markersize=legend.mark_size, label=str(val)))
            labels.append(str(val))

    # COLORS #
    if style.hue_col:
        current_palette = sns.color_palette()
        hue_vals = sorted(df_data[style.hue_col].unique())
        for idx, val in enumerate(hue_vals):
            col = style.color_map.get(val, current_palette[idx]) if style.color_map else current_palette[idx]
            handles.append(plt.Line2D([0], [0], marker="o", linestyle="", color=col, markersize= legend.mark_size, label=str(val)))
            labels.append(str(val))

    # ----- LEGEND OUTPUT AND PLACEMENT ----- #

    if handles:
        if legend.path:
            leg_fig, leg_ax = plt.subplots(figsize=(2, 2))
            leg_ax.axis("off")
            leg_ax.legend(handles=handles, labels=labels, loc="center", ncol=legend.ncol, frameon=legend.frameon, title=legend.title)
            Path(legend.path).parent.mkdir(parents=True, exist_ok=True)
            path = Path(f"{str(legend.path)}.{output.out_type}")
            leg_fig.savefig(path, dpi=output.dpi, transparent=output.transparent, format=output.out_type)
            plt.close(leg_fig)

        if legend.loc:
            fig.legend(handles=handles, labels=labels, loc=legend.loc,
                   ncol=legend.ncol, frameon=legend.frameon,
                   title=legend.title)

    # ----- OUTPUT ----- #

    # RASTERIZE #
    if output.rasterize:
        for ax in axes:
            for coll in ax.collections:
                coll.set_rasterized(True)
    plt.tight_layout()
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

    return fig, axes
