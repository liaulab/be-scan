
from be_scan.figure_plot.figure_classes import *

def _style_attr(style, name):
    """
    Read a style field, falling back to the PWES defaults.

    pwes_histogram and pwes_scatterplots are annotated
    `style: ScatterStyleOpts = SCATTERSTYLE` but read linewidth, edgecolor and
    the spine_* fields, which only exist on PWESHeatmapStyleOpts. With their
    own default they raised AttributeError, so neither could be called without
    passing a style object of the other type. Reading through here lets either
    dataclass work.
    """
    return getattr(style, name, getattr(PWESHEATMAPSTYLE, name))


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
        linewidth=_style_attr(style, 'linewidth'), edgecolor=_style_attr(style, 'edgecolor'))

    plt.xlim(min(clust_counts.index)-1, max(clust_counts.index)+1)
    plt.xticks(clust_counts.index, fontproperties=arial_font6)
    plt.ylim(axis.ylim)
    plt.yticks(fontproperties=arial_font6)

    plt.xlabel("Cluster Number", fontproperties=arial_font6)
    plt.ylabel("Count", fontproperties=arial_font6)
    plt.title("Number of Guides in Clusters", fontproperties=arial_font6)

    finish_figure(plt.gcf(), output, suffix="-cluster_histogram")

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
        **_style_attr(style, 'swarmplot_kws'),
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

    finish_figure(plt.gcf(), output, suffix="-cluster_swarm", bbox_inches="tight")

    # SWARMPLOT BOXPLOT #

    fig, ax = plt.subplots(figsize=output.figsize)
    ax = plt.gca()

    sns.swarmplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        palette=color_list, zorder=1,
        **_style_attr(style, 'swarmplot_kws'),
        )
    for collection in ax.collections:
        collection.set_rasterized(True)
    sns.boxplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        zorder=2,
        **_style_attr(style, 'boxplot_kws'),
        boxprops=_style_attr(style, 'boxprops'),
        whiskerprops=_style_attr(style, 'whiskerprops'),
        capprops=_style_attr(style, 'capprops'),
        medianprops=_style_attr(style, 'medianprops'),
        )
    ax.tick_params(width=0.5)

    plt.xlim(min(df_clus['cl_new'])-2, max(df_clus['cl_new']))
    plt.ylim(axis.ylim)
    plt.xticks(fontproperties=arial_font6)
    plt.yticks(fontproperties=arial_font6)

    plt.ylabel(f'{x_col} Values', fontproperties=arial_font6)
    plt.xlabel('Cluster Number', fontproperties=arial_font6)
    plt.title('Clusters Boxplot', fontproperties=arial_font6)

    finish_figure(plt.gcf(), output, suffix="-cluster_swarm_boxplot", bbox_inches="tight")

    # STRIPPLOT #

    fig, ax = plt.subplots(figsize=output.figsize)
    ax = plt.gca()

    sns.stripplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        palette=color_list, zorder=1,
        jitter=True,
        **_style_attr(style, 'stripplot_kws'),
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

    finish_figure(plt.gcf(), output, suffix="-cluster_strip", bbox_inches="tight")

    # STRIPPLOT BOXPLOT #

    fig, ax = plt.subplots(figsize=output.figsize)
    ax = plt.gca()

    sns.stripplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        palette=color_list, zorder=1,
        jitter=True,
        **_style_attr(style, 'stripplot_kws'),
        )
    for collection in ax.collections:
        collection.set_rasterized(True)
    sns.boxplot(
        ax=ax,
        x='cl_new', y=x_col, data=df_clus,
        zorder=2,
        **_style_attr(style, 'boxplot_kws'),
        boxprops=_style_attr(style, 'boxprops'),
        whiskerprops=_style_attr(style, 'whiskerprops'),
        capprops=_style_attr(style, 'capprops'),
        medianprops=_style_attr(style, 'medianprops'),
        )
    ax.tick_params(width=0.5)

    plt.xlim(min(df_clus['cl_new'])-2, max(df_clus['cl_new']))
    plt.ylim(axis.ylim)
    plt.xticks(fontproperties=arial_font6)
    plt.yticks(fontproperties=arial_font6)

    plt.ylabel(f'{x_col} Values', fontproperties=arial_font6)
    plt.xlabel('Cluster Number', fontproperties=arial_font6)
    plt.title('Clusters Boxplot', fontproperties=arial_font6)

    finish_figure(plt.gcf(), output, suffix="-cluster_strip_boxplot", bbox_inches="tight")

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
        for dom in resolve_domains(domain):
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
            # wrap, so fewer colours than clusters cycles rather than raising
            color=color_list[(cluster-1) % len(color_list)],
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
        ax.set_title(f'Cluster {cluster}', fontproperties=arial_font6)

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
        **_style_attr(style, 'heatmap_kws'),
        )

    for edge, spine in ax.spines.items():
        spine.set_visible(_style_attr(style, 'spine_visible'))
        spine.set_color(_style_attr(style, 'spine_color'))

    # SET TICKS AND LABELS #
    ax.set_xlabel(axis.xlabel, fontproperties=arial_font6)
    ax.set_ylabel(axis.ylabel, fontproperties=arial_font6)
    if axis.xticks: ax.set_xticklabels(axis.xticks)
    if axis.yticks: ax.set_yticklabels(axis.yticks)
    # TICK LABELS #
    for label in ax.get_xticklabels() + ax.get_yticklabels(): label.set_fontproperties(arial_font6)
    plt.subplots_adjust(wspace=_style_attr(style, 'spine_wspace'), hspace=_style_attr(style, 'spine_hspace'))

    # DRAW RECTANGLES #
    # HORIZONTAL #
    for i, idx in enumerate(df_scaled.index):
        for start in blackout:
            # labels are integer residue numbers for a single subunit and
            # gene-prefixed strings for a complex; only the latter has a prefix
            if str(idx).startswith(str(start)):
                rect = patches.Rectangle(
                    (0, i), len(df_scaled.columns), 1,
                    linewidth=0, edgecolor=None, facecolor='black', alpha=1)
                ax.add_patch(rect)
    # VERTICAL #
    for j, col in enumerate(df_scaled.columns):
        for start in blackout:
            if str(col).startswith(str(start)):
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
        df_colors['Colors'] = [color_list[(i-1) % len(color_list)]
                               for i in df_clusters['cl_new']]
        rcolor = df_colors['Colors'].copy()
    else: rcolor = None

    # CLUSTERMAP #
    fig = sns.clustermap(
        df_scaled,
        row_linkage=link, col_linkage=link,
        row_colors=rcolor,
        figsize=output.figsize,
        **_style_attr(style, 'clustered_heatmap_kws'),
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
        spine.set_visible(_style_attr(style, 'spine_visible'))
        spine.set_color(_style_attr(style, 'spine_color'))

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

    norm = mpl.colors.Normalize(**_style_attr(style, 'colorbar_normalize_kws'))
    cbar = mpl.colorbar.ColorbarBase(
        ax, norm=norm, orientation=orientation,
        **_style_attr(style, 'colorbar_base_kws'))

    for label in (ax.get_yticklabels() if orientation == 'vertical' else ax.get_xticklabels()):
        label.set_fontproperties(arial_font6)

    ax.tick_params(**axis.tick_kws)
    cbar.outline.set_linewidth(_style_attr(style, 'colorbar_linewidth'))

    if output.save:
        out_path = Path(f"{str(output.path)}.{output.out_type}")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, format=output.out_type,
                    dpi=output.dpi, transparent=output.transparent
                    )
    # SHOW #
    if output.show: plt.show()
    plt.close(fig)

# PWES Cluster Mapping #

def tuple_to_hex(rgb_tuple):
    if not (isinstance(rgb_tuple, tuple) and len(rgb_tuple) == 3 and
            all(isinstance(val, float) and 0 <= val <= 1 for val in rgb_tuple)):
        return None

    r, g, b = rgb_tuple
    r, g, b = int(r * 255), int(g * 255), int(b * 255)
    return f"{r:02x}{g:02x}{b:02x}".upper()

def generate_pymol_script(
    output_filename,
    pdb_filename,
    residue_dict,
    colors, 
    ):

    """
    Generates a PyMOL script to display selected residues as separate objects.
    """

    with open(output_filename, 'w') as pymol_script:
        pymol_script.write(f"load {pdb_filename}\n")
        pymol_script.write(f"show cartoon\n")
        pymol_script.write(f"color white\n")
        pymol_script.write(f"zoom\n")

        for i, (key, vals) in enumerate(residue_dict.items()):
            group_name = f"Cluster{key}"
            # res has a format of ('A', 123) #
            selection_str = " or ".join([f"resi {str(int(res[1:]))} and chain {res[0]} and name CA" for res in vals])

            color = colors[i]
            color_hex = tuple_to_hex(color)
            color_str = f'0x{color_hex}'
            print(color_str)
            
            pymol_script.write(f"select {selection_str}\n")
            pymol_script.write(f"create {group_name}, sele\n")
            pymol_script.write(f"show spheres, {group_name}\n")
            pymol_script.write(f"color {color_str}, {group_name}\n")
            pymol_script.write(f"disable {group_name}\n")
