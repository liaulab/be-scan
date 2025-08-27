
from be_scan.figure_plot.figure_classes import *

def pwes_histogram(
    df_clus,
    color_list,
    
    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    """
    Given df_clus, generated from cluster_pws, plot histogram of the number of guides in each cluster.
    """

    assert "cl_new" in df_clus.columns, "df_clus does not contain a 'cl_new' column"
    clust_counts = df_clus["cl_new"].value_counts().sort_index()

    # HISTOGRAM #

    fig, ax = plt.subplots(figsize=output.figsize)
    plt.bar(
        clust_counts.index, clust_counts,
        color=color_list,
        linewidth=style.linewidth, edgecolor=style.edgecolor)

    plt.xlim(min(clust_counts.index)-1, max(clust_counts.index)+1)
    plt.xticks(clust_counts.index, fontproperties=arial_font6)
    plt.ylim(axis.ylim)
    plt.yticks(fontproperties=arial_font6)

    plt.xlabel("Cluster Number", fontproperties=arial_font6)
    plt.ylabel("Count", fontproperties=arial_font6)
    plt.title("Number of Guides in Clusters", fontproperties=arial_font6)

    out_path = Path(f"{str(output.path)}-cluster_histogram.{output.out_type}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(
        out_path, dpi=output.dpi,
        format=output.out_type, transparent=output.transparent,
        )
    if output.show: plt.show()
    plt.close()

def pwes_boxplots(
    df_clus,
    x_col,
    color_list,
    
    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    """
    Given df_clus, generated from cluster_pws,
    plot boxplot variations of the values for each cluster of guides.
    """

    assert "cl_new" in df_clus.columns, "df_clus does not contain a 'cl_new' column"
    assert x_col in df_clus.columns, f"{x_col} not in df_clus"

    # SWARMPLOT #

    fig, ax = plt.subplots(figsize=output.figsize)
    ax = plt.gca()

    sns.swarmplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        palette=color_list, zorder=1,
        **style.swarmplot_kws,
        )
    for collection in ax.collections:
        collection.set_rasterized(True)
    ax.tick_params(width=0.5)

    plt.xlim(min(df_clus['cl_new'])-2, max(df_clus['cl_new']))
    plt.ylim(axis.ylim)
    plt.xticks(fontproperties=arial_font6)
    plt.yticks(fontproperties=arial_font6)

    plt.ylabel(f'{x_col} Values', fontproperties=arial_font6)
    plt.xlabel('Cluster Number', fontproperties=arial_font6)
    plt.title('Clusters Boxplot', fontproperties=arial_font6)

    out_path = Path(f"{str(output.path)}-cluster_swarm.{output.out_type}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=output.dpi,
                format=output.out_type, transparent=output.transparent,
                bbox_inches="tight")
    if output.show: plt.show()
    plt.close()

    # SWARMPLOT BOXPLOT #

    fig, ax = plt.subplots(figsize=output.figsize)
    ax = plt.gca()

    sns.swarmplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        palette=color_list, zorder=1,
        **style.swarmplot_kws,
        )
    for collection in ax.collections:
        collection.set_rasterized(True)
    sns.boxplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        zorder=2,
        **style.boxplot_kws,
        boxprops=style.boxprops,
        whiskerprops=style.whiskerprops,
        capprops=style.capprops,
        medianprops=style.medianprops,
        )
    ax.tick_params(width=0.5)

    plt.xlim(min(df_clus['cl_new'])-2, max(df_clus['cl_new']))
    plt.ylim(axis.ylim)
    plt.xticks(fontproperties=arial_font6)
    plt.yticks(fontproperties=arial_font6)

    plt.ylabel(f'{x_col} Values', fontproperties=arial_font6)
    plt.xlabel('Cluster Number', fontproperties=arial_font6)
    plt.title('Clusters Boxplot', fontproperties=arial_font6)

    out_path = Path(f"{str(output.path)}-cluster_swarm_boxplot.{output.out_type}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=output.dpi,
                format=output.out_type, transparent=output.transparent,
                bbox_inches="tight")
    if output.show: plt.show()
    plt.close()

    # STRIPPLOT #

    fig, ax = plt.subplots(figsize=output.figsize)
    ax = plt.gca()

    sns.stripplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        palette=color_list, zorder=1,
        jitter=True,
        **style.stripplot_kws,
        )
    for collection in ax.collections:
        collection.set_rasterized(True)
    ax.tick_params(width=0.5)

    plt.xlim(min(df_clus['cl_new'])-2, max(df_clus['cl_new']))
    plt.ylim(axis.ylim)
    plt.xticks(fontproperties=arial_font6)
    plt.yticks(fontproperties=arial_font6)

    plt.ylabel(f'{x_col} Values', fontproperties=arial_font6)
    plt.xlabel('Cluster Number', fontproperties=arial_font6)
    plt.title('Clusters Boxplot', fontproperties=arial_font6)

    out_path = Path(f"{str(output.path)}-cluster_strip.{output.out_type}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=output.dpi,
                format=output.out_type, transparent=output.transparent,
                bbox_inches="tight")
    if output.show: plt.show()
    plt.close()

    # STRIPPLOT BOXPLOT #

    fig, ax = plt.subplots(figsize=output.figsize)
    ax = plt.gca()

    sns.stripplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        palette=color_list, zorder=1,
        jitter=True,
        **style.stripplot_kws,
        )
    for collection in ax.collections:
        collection.set_rasterized(True)
    sns.boxplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        zorder=2,
        **style.boxplot_kws,
        boxprops=style.boxprops,
        whiskerprops=style.whiskerprops,
        capprops=style.capprops,
        medianprops=style.medianprops,
        )
    ax.tick_params(width=0.5)

    plt.xlim(min(df_clus['cl_new'])-2, max(df_clus['cl_new']))
    plt.ylim(axis.ylim)
    plt.xticks(fontproperties=arial_font6)
    plt.yticks(fontproperties=arial_font6)

    plt.ylabel(f'{x_col} Values', fontproperties=arial_font6)
    plt.xlabel('Cluster Number', fontproperties=arial_font6)
    plt.title('Clusters Boxplot', fontproperties=arial_font6)

    out_path = Path(f"{str(output.path)}-cluster_strip_boxplot.{output.out_type}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=output.dpi,
                format=output.out_type, transparent=output.transparent,
                bbox_inches="tight")
    if output.show: plt.show()
    plt.close()

def pwes_scatterplots(
    df_clus,
    x_col, scores_col,
    color_list,
    max_cluster=None,

    axis: AxisLabelOpts = AXIS,
    style: ScatterStyleOpts = SCATTERSTYLE,
    domain: DomainOpts = DOMAIN,
    output: OutputOpts = OUTPUT,
    ):

    """
    Generates separate scatter plots for each cluster in `cl_new`,
    highlighting one cluster per plot while coloring others gray.
    """

    try: df_clus['cl_new']
    except KeyError: raise Exception('df_clus does not contain a "cl_new" column')

    unique_clusters = sorted(df_clus['cl_new'].unique())
    if max_cluster is None: max_cluster = max(unique_clusters)

    for cluster in range(max_cluster):
        cluster = int(cluster+1)

        fig, ax = plt.subplots(1, 1, figsize=output.figsize) # SUBPLOTS #

        # PLOT DOMAINS #
        for dom in domain.dom_setting:
            ax.axvspan(
                dom["start"], dom["end"],
                facecolor=dom.get("color", domain.default_color),
                edgecolor="none",
                alpha=dom.get("alpha", domain.default_alpha) )

        # SCATTERPLOT GRAY FOR ALL OTHER POINTS #
        sns.scatterplot(
            data=df_clus[df_clus['cl_new'] != cluster],
            ax=ax, legend=False,
            x=x_col, y=scores_col,
            **style.transparent_kws,
        )
        # SCATTERPLOT HIGHLIGHT CURRENT CLUSTER #
        sns.scatterplot(
            data=df_clus[df_clus['cl_new'] == cluster],
            ax=ax,
            x=x_col, y=scores_col,
            color=color_list[cluster-1],
            **style.opaque_kws,
        )

        # SET LIMITS #
        if axis.ylim: ax.set_ylim(*axis.ylim)
        if axis.xlim: ax.set_xlim(*axis.xlim)

        ax.spines['left'].set_position(('outward', 3))
        ax.spines['right'].set_position(('outward', 3))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        fig_leg = ax.get_legend()
        if fig_leg is not None: fig_leg.remove()

        # SET X-TICKS #
        if axis.xticks: ax.set_xticks(axis.xticks)
        elif axis.protein_len: ax.set_xticks([1, axis.protein_len])
        # SET Y-TICKS #
        if axis.yticks: ax.set_yticks(axis.yticks)
        # TICK LABELS #
        for label in ax.get_xticklabels() + ax.get_yticklabels(): label.set_fontproperties(arial_font6)

        # AXES LABELS #
        ax.set_xlabel(None, fontproperties=arial_font6)
        ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)

        # ----- OUTPUT ----- #

        # RASTERIZE #
        if output.rasterize and output.rasterize == True:
            for coll in ax.collections:
                coll.set_rasterized(True)

        plt.tight_layout()
        if output.save:
            out_path = Path(f"{str(output.path)}_subscatter_Cluster{str(cluster)}.{output.out_type}")
            out_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_path, format=output.out_type,
                        dpi=output.dpi, transparent=output.transparent
                        )
        # SHOW #
        if output.show: plt.show()
        plt.close(fig)

def pwes_heatmap(
    df_scaled,
    mask_on=True,
    blackout=['X'],

    axis: AxisLabelOpts = AXIS,
    style: PWESHeatmapStyleOpts = PWESHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    # SETUP #
    mask = np.zeros_like(df_scaled)
    if mask_on: mask[np.tril_indices_from(mask, k=-1)] = True # MASK FOR LOWER TRIANGLE #

    # GENERATE HEATMAP #
    fig, ax = plt.subplots(figsize=output.figsize)
    sns.heatmap(
        data=df_scaled,
        ax=ax,
        mask=mask,
        **style.heatmap_kws,
        )

    for edge, spine in ax.spines.items():
        spine.set_visible(style.spine_visible)
        spine.set_color(style.spine_color)

    # SET TICKS AND LABELS #
    ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    if axis.xticks: ax.set_xticklabels(axis.xticks)
    if axis.yticks: ax.set_yticklabels(axis.yticks)
    # TICK LABELS #
    for label in ax.get_xticklabels() + ax.get_yticklabels(): label.set_fontproperties(arial_font6)
    plt.subplots_adjust(wspace=style.spine_wspace, hspace=style.spine_hspace)

    # DRAW RECTANGLES #
    # HORIZONTAL #
    for i, idx in enumerate(df_scaled.index):
        for start in blackout:
            if idx.startswith(start):
                rect = patches.Rectangle(
                    (0, i), len(df_scaled.columns), 1,
                    linewidth=0, edgecolor=None, facecolor='black', alpha=1)
                ax.add_patch(rect)
    # VERTICAL #
    for j, col in enumerate(df_scaled.columns):
        for start in blackout:
            if col.startswith(start):
                rect = patches.Rectangle(
                    (j, 0), 1, len(df_scaled),
                    linewidth=0, edgecolor=None, facecolor='black', alpha=1)
                ax.add_patch(rect)

    # RASTERIZE #
    if output.rasterize:
        for coll in ax.collections:
            coll.set_rasterized(True)

    plt.tight_layout()
    if output.save:
        out_path = Path(f"{str(output.path)}_heatmap.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent
                    )
    # SHOW #
    if output.show: plt.show()
    plt.close(fig)

def pwes_clustermap(
    df_scaled, link, df_clusters,
    color_list=[],

    axis: AxisLabelOpts = AXIS,
    style: PWESHeatmapStyleOpts = PWESHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    # SET COLORS FOR CLUSTERS #
    if color_list:
        df_colors = pd.DataFrame(index=df_scaled.index, columns=['Cluster'])
        df_colors['Colors'] = [color_list[i-1] for i in df_clusters['cl_new']]
        rcolor = df_colors['Colors'].copy()
    else: rcolor = None

    # CLUSTERMAP #
    fig = sns.clustermap(
        df_scaled,
        row_linkage=link, col_linkage=link,
        row_colors=rcolor,
        figsize=output.figsize,
        **style.clustered_heatmap_kws,
        )
    ax = fig.ax_heatmap

    # AXIS SETTINGS #
    fig.ax_heatmap.collections[0].set_rasterized(output.rasterize)
    fig.cax.set_visible(False)
    ax.set_aspect('equal')
    fig.ax_col_dendrogram.set_visible(False)

    # DRAW LINES BETWEEN CLUSTERS #
    temp = 0
    clusterlist = df_clusters['cl_new'].value_counts().sort_index()
    for i in clusterlist.iloc[:-1]:
        temp += i
        ax.axhline(y=temp, **style.line_kws)

    for edge, spine in ax.spines.items():
        spine.set_visible(style.spine_visible)
        spine.set_color(style.spine_color)

    plt.tight_layout()

    # SAVE #
    if output.save:
        out_path = Path(f"{str(output.path)}_cluster_heatmap.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent
                    )
    # SHOW #
    if output.show: plt.show()
    plt.close()

def pwes_heatmap_colorbar(
    orientation='vertical', # 'v' or 'h'

    axis: AxisLabelOpts = AXIS,
    style: PWESHeatmapStyleOpts = PWESHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    fig, ax = plt.subplots(figsize=output.figsize)

    norm = mpl.colors.Normalize(**style.colorbar_normalize_kws)
    cbar = mpl.colorbar.ColorbarBase(
        ax, norm=norm, orientation=orientation,
        **style.colorbar_base_kws)

    for label in (ax.get_yticklabels() if orientation == 'vertical' else ax.get_xticklabels()):
        label.set_fontproperties(arial_font6)

    ax.tick_params(**axis.tick_kws)
    cbar.outline.set_linewidth(style.colorbar_linewidth)

    if output.save:
        out_path = Path(f"{str(output.path)}.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent
                    )
    # SHOW #
    if output.show: plt.show()
    plt.close(fig)
