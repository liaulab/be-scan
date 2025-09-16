
from be_scan.figure_plot.figure_classes import *

def correlation_scatterplot_figure(
    df_data, # dataframe
    x_column, # the x axis values
    y_column, # the y axis values
    legend_order,

    hue_col='',
    palette = {},
    highlight_col="mut_type",
    highlight_categories=["Splice", "Intron", "Nonsense"],

    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    """[Summary]
    This function takes in a dataframe from analysis and then
    plots the data for each condition to reveal which guides are enriched
    """
    fig, ax = plt.subplots(1, 1, figsize=output.figsize)

    # OPACITY COLUMN #
    thr = style.transparency_threshold
    upper, lower = -thr, thr
    df_data['_alpha'] = np.where(
        (df_data[x_column].between(upper, lower)) & (df_data[y_column].between(upper, lower)),
        "low", "high",
    )

    # CUTOFF LINES #
    for cutoff in [-thr, thr]:
        ax.axhline(y=cutoff, **style.cutoff_line_kws)
        ax.axvline(x=cutoff, **style.cutoff_line_kws)

    # OUTLINE SPLICE AND INTRONS #
    df_data_si = df_data.loc[df_data[highlight_col].isin(highlight_categories)]
    df_data_no_si = df_data.loc[~df_data[highlight_col].isin(highlight_categories)]

    # SCATTER NOT HIGHLIGHTED CATEGORY #
    df_data_temp = df_data_no_si[df_data_no_si['_alpha'] == "high"]
    scatter = sns.scatterplot(
        ax=ax, data=df_data_temp, x=x_column, y=y_column,
        hue=hue_col, palette = palette,
        **style.opaque_kws, ) # alpha=1.0, linewidth=0, s=4
    df_data_temp = df_data_no_si[df_data_no_si['_alpha'] == "low"]
    scatter = sns.scatterplot(
        ax=ax, data=df_data_temp, x=x_column, y=y_column,
        hue=hue_col, palette = palette,
        **style.transparent_kws, ) # alpha=0.25, linewidth=0, s=4

    # SCATTER OUTLINE HIGHLIGHTED CATEGORY #
    df_data_temp = df_data_si[df_data_si['_alpha'] == "high"]
    scatter = sns.scatterplot(
        ax=ax, data=df_data_temp, x=x_column, y=y_column,
        hue=hue_col, palette = palette,
        **style.highlight_kws) #alpha=1.0, linewidth=0.3, edgecolor='black', s=4

    if output.rasterize:
        for coll in scatter.collections: coll.set_rasterized(True)

    # LINE OF BEST FIT
    x_vals = df_data[x_column].values
    y_vals = df_data[y_column].values

    mask = ~np.isnan(x_vals) & ~np.isnan(y_vals)
    slope, intercept = np.polyfit(x_vals[mask], y_vals[mask], 1)
    print(f'Best fit line: y = {slope:.4f}x + {intercept:.4f}')

    x_fit = np.linspace(np.min(x_vals), np.max(x_vals), 100)
    y_fit = slope * x_fit + intercept
    ax.plot(x_fit, y_fit, **style.line_kws)

    # PEARSON #
    r_val, p_val = pearsonr(x_vals[mask], y_vals[mask])
    print(f'Pearson r = {r_val:.4f}, p = {p_val:.4e}')

    # AXES LABELS AND TICKS #
    if axis.title: ax.set_title(axis.title, fontproperties=arial_font6)
    if axis.xlabel: ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    if axis.ylabel: ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    if axis.xlim: ax.set_xlim(*axis.xlim)
    if axis.ylim: ax.set_ylim(*axis.ylim)
    if axis.xticks: ax.set_xticks(axis.xticks)
    if axis.yticks: ax.set_yticks(axis.yticks)
    for label in ax.get_xticklabels() + ax.get_yticklabels(): label.set_fontproperties(arial_font6)

    ax.spines['left'].set_position(('outward', 3))
    ax.spines['left'].set_linewidth(axis.linewidth)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig_leg = ax.get_legend()
    if fig_leg is not None: fig_leg.remove()

    # LEGEND #
    available_labels = df_data[hue_col].dropna().unique().tolist()
    final_order = [label for label in legend_order if label in available_labels]
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='k', label=label,
                   markerfacecolor=palette[label], markersize=4, linewidth=0)
        for label in final_order if label in palette
    ]

    # PLOT LEGEND #
    if output.save:
        fig_legend = plt.figure(figsize=output.figsize)
        fig_legend.legend(handles=legend_elements, loc='center', frameon=False)
        fig_legend.tight_layout()
        legend_out = Path(f"{str(output.path)}-legend.{output.out_type}")
        legend_out.parent.mkdir(parents=True, exist_ok=True)
        fig_legend.savefig(legend_out, format=output.out_type,
                           dpi=output.dpi, transparent=output.transparent)
    plt.close(fig_legend)

    # SAVE #
    if output.save:
        out_path = Path(f"{str(output.path)}.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent)
    if output.show: plt.show()
    plt.close(fig)
