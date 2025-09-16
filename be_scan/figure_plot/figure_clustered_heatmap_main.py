
from be_scan.figure_plot.figure_classes import *

def kmeans_clustered_heatmap(
    df, k,
    clusters, # CLUSTERS CORRESPONDING TO EACH GUIDE
    genes,
    screen_order=None,
    gene_palette={},
    gene_list=None, # LIST OF GENES PRESENT IN SCREEN
    cluster_order=None,
    highlight_cols=None,
    order_col=None,

    axis: AxisLabelOpts = AXIS,
    style: ClusteredHeatmapStyleOpts = CLUSTEREDHEATMAPSTYLE,
    output: OutputOpts = OUTPUT,
    ):

    # REORDER BY SCREENS (COLUMNS)
    if screen_order:
        for screen in screen_order:
            assert screen in df.columns, f'{screen} not in dataframe'
        df = df[screen_order]
    else: screen_order = df.columns.tolist()

    # CLUSTER ASSIGNMENTS FROM PREVIOUS CLUSTERING AND ORDER #
    df['Clusters'] = clusters
    df['Clusters'] = df['Clusters'].astype("category")
    df['Genes'] = genes
    sorted_df = df.sort_values("Clusters")

    # ORDER BY ORDER COL IF PROVIDED #
    if order_col is not None and order_col in df.columns:
        # SORT BY CLUSTER THEN VALUE #
        sorted_df["Clusters"] = pd.Categorical(sorted_df["Clusters"], categories=cluster_order, ordered=True)
        sorted_df = sorted_df.sort_values(["Clusters", order_col], ascending=[True, False])
    elif order_col is not None:
        print('[order_col] not found in dataframe.')
    else:
        # SORT BY CLUSTER ONLY #
        sorted_df = sorted_df.sort_values("Clusters")
    sorted_rows = sorted_df.index.tolist()

    # COLOR #
    if not style.palette: palette = sns.color_palette(style.color_pal, n_colors=k)
    else: palette = style.palette
    sorted_df['Clusters'] = sorted_df['Clusters'].astype(int)
    sorted_df['Colors'] = sorted_df['Clusters'].map(lambda x: palette[x])
    cluster_colors_array = np.array(sorted_df['Colors'].tolist()).reshape(-1, 1, 3)

    # GENES, HIGHLIGHT, CLUSTERS #
    total_cols = 1
    if gene_list:
        total_cols += len(gene_list)
    if highlight_cols:
        total_cols += len(highlight_cols)

    # CREATE FIGURE #
    fig = plt.figure(
        figsize=output.figsize,
        )
    gs = fig.add_gridspec(
        nrows=1, ncols=total_cols+1,
        width_ratios=[style.col_width]*total_cols + [style.heatmap_width],
        wspace=style.wspace,
        )

    i = 0
    # GENE COLOR BARS #
    gene_palette['Default'] = style.default_color
    if gene_list:
        for gene in gene_list:
            gene_assignments = pd.Series(
                [g if gene == g else 'Default' for g in sorted_df['Genes']],
                index=sorted_df.index,
                )
            sorted_genes = gene_assignments.loc[sorted_rows] # SORT ORDER #
            gene_colors = sorted_genes.map(
                lambda x: gene_palette.get(x, (1, 1, 1))
                )
            gene_color_array = np.array(gene_colors.tolist()).reshape(-1, 1, 3)
            ax_gene = fig.add_subplot(gs[0, i])
            ax_gene.imshow(gene_color_array, aspect="auto")
            ax_gene.axis("off")
            i += 1
    # HIGHLIGHT BARS #
    if highlight_cols:
        for highlight_col in highlight_cols:
            highlight_assignments = sorted_df['Genes'].where(
                sorted_df.index.isin(highlight_col), 'Default'
            )
            sorted_highlights = highlight_assignments.loc[sorted_rows] # SORT ORDER #
            highlight_colors = sorted_highlights.map(
                lambda x: gene_palette.get(x, (1, 1, 1))
                )
            highlight_color_arr = np.array(highlight_colors.tolist()).reshape(-1, 1, 3)
            ax_highlight = fig.add_subplot(gs[0, i])
            ax_highlight.imshow(highlight_color_arr, aspect="auto")
            ax_highlight.axis("off")
            i += 1
    # CLUSTER BAR #
    ax_cluster = fig.add_subplot(gs[0, i])
    ax_cluster.imshow(cluster_colors_array, aspect="auto")
    ax_cluster.axis("off")
    i += 1
    # HEATMAP #
    ax_heatmap = fig.add_subplot(gs[0, i])
    sns.heatmap(
        sorted_df[screen_order],
        ax=ax_heatmap,
        **style.heatmap_kws,
        )

    # ADD BORDER #
    num_cols, num_rows = len(screen_order), sorted_df.shape[0]
    ax_heatmap.add_patch(Rectangle(
        (0, 0), num_cols, num_rows,
        **style.border_kws,
    ))
    for i in range(1, num_cols):
        ax_heatmap.axvline(i, **style.line_kws)

    # LABELS #
    ax_heatmap.set_title(axis.title, fontproperties=arial_font6)
    ax_heatmap.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    ax_heatmap.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    ax_heatmap.set_yticks(axis.yticks)
    for label in ax_heatmap.get_xticklabels():
        label.set_fontproperties(arial_font6)

    # SAVE #
    out_path = Path(f"{str(output.path)}.{output.out_type}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=output.dpi,
                format=output.out_type, transparent=output.transparent)
    plt.show()
    plt.close()
