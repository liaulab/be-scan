
from be_scan.figure_plot.figure_classes import *
from be_scan.figure_plot.figure_pwes_helpers import *

def pwes_clustering(
    df_scores,
    x_col, scores_col, pdb_file,
    gene_col='', gene_map={},
    tanh_a=1, gauss_std = 16, dend_t = 14,
    aa_int=None, out_prefix=None, out_dir=None,
    pos_only=False,
    ):

    """
    Main function to run 3D clustering analysis.
    """
    # BASIC CHECKS #
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    assert x_col in df_scores.columns, 'Check [x_col] in df_scores'
    assert scores_col in df_scores.columns, 'Check [scores_col] in df_scores'

    # CALCULATE SPATIAL COMPONENT FROM STRUCTURE #
    df_centroids = process_pdb(pdb_file, gene_map.values())
    df_pwdist = get_pairwise_dist(df_centroids, gene_map.values(), aa_int)
    df_gauss = df_pwdist.apply(lambda x: gauss(x, gauss_std))

    # ASSIGN IDs FOR EVERY INPUT #
    sgrnaID = ["sgRNA_" + num for num in map(str, list(range(df_scores.shape[0])))]
    df_scores["sgRNA_ID"] = sgrnaID

    df_scores[x_col] = df_scores[x_col].astype(int) # MUST BE INTEGER #
    # FOR ONLY ONE SUBUNIT #
    if len(gene_map) == 0 or gene_col == 0:
        print('No chain(s) indicated. Automatically assigning to structure ...')
        aa_in_structure = df_centroids.aa_num.tolist()
        df_scores = df_scores[df_scores[x_col].isin(aa_in_structure)]

        # ONLY USE BASE EDITING SCORES FOR RESIDUES FOUND IN COMPLEX #
        list_aas = df_scores[x_col]
        df_pws_score = calculate_pw_score(df_scores, scores_col, tanh_a)
    # FOR A COMPLEX #
    else:
        print('Chain(s) indicated. Mapping gene(s) to chains ...')
        assert gene_col in df_scores.columns, 'Check [gene_col] is in '

        aa_in_structure = df_centroids.label.tolist()
        df_scores = df_scores[df_scores[gene_col].isin(gene_map.keys())] # FILTER OUT GENES NOT IN MAP #
        df_scores['label'] = df_scores[gene_col].map(gene_map) + df_scores[x_col].astype(str).str.zfill(4) # FORMAT GENE-POSITION
        df_scores = df_scores[df_scores['label'].isin(aa_in_structure)]

        list_aas = df_scores['label']
        df_pws_score = calculate_pw_score(df_scores, scores_col, tanh_a)

    del df_centroids, df_pwdist

    # CALCULATE PWES SCORE #
    df_pwes_sorted, df_pwes_unsorted = calculate_pwes(df_gauss, df_pws_score, list_aas, pos_only)

    # CLUSTER PWES #
    if len(gene_map) == 0 or gene_col == 0: cluster_xcol = x_col
    else: cluster_xcol = 'label'
    df_clus, link = cluster_pws(
        df_pws = df_pwes_unsorted, df_score = df_scores,
        list_aas=list_aas, t=dend_t, x_col=cluster_xcol)

    # SAVE CLUSTERS #
    aas_dict = get_clus_aa(df_clus, cluster_xcol)
    with open(f"{out_dir}{out_prefix}_aas_dict.pickle", "wb") as file:
        pickle.dump(aas_dict, file)

    # SAVE DATAFRAMES #
    out = f"{out_dir}{out_prefix}"
    Path(f"{out}").parent.mkdir(parents=True, exist_ok=True)
    df_gauss.to_csv(f"{out}_df_gauss.csv")
    df_pwes_sorted.to_csv(f"{out}_df_pwes_sorted.csv")
    df_pwes_unsorted.to_csv(f"{out}_df_pwes_unsorted.csv")
    df_clus.to_csv(f"{out}_df_clus.csv")

    return df_gauss, df_pwes_sorted, df_pwes_unsorted, df_clus, aas_dict, link
