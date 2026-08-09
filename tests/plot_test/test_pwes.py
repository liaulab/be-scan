"""
PWES 3D clustering.

Rewritten for 3.0.0. The previous version read its CSVs at import time, so a
missing file broke collection of the entire suite rather than skipping these
tests, and it called be_scan.plot.pwes_clustering(show_plots=False) -- a
module and a parameter that do not exist. pwes_clustering lives in
be_scan.figure_plot and takes no show_plots.

Data is loaded inside fixtures and skipped cleanly when absent, so the suite
still runs on a checkout without the .pdb files.
"""

import pickle

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import pytest

from be_scan.figure_plot import pwes_clustering


DATA = "tests/test_data/plot"
SCORES = f"{DATA}/q575r_ezh2_abe_missense_pos.csv"
AF_PDB = f"{DATA}/AF_Q15910.pdb"
XTAL_PDB = f"{DATA}/6wkr-chainC-EZH2.pdb"

PWES_PARAMS = dict(x_col="xpos", scores_col="DMSO-PRC2i",
                   gauss_std=16, dend_t=14, tanh_a=1.0)


@pytest.fixture(scope="module")
def scores(request):
    path = request.config.rootpath / SCORES
    if not path.exists():
        pytest.skip(f"{SCORES} not available")
    return pd.read_csv(path)


@pytest.fixture(params=["alphafold", "crystal"])
def pdb(request):
    path = request.config.rootpath / (
        AF_PDB if request.param == "alphafold" else XTAL_PDB)
    if not path.exists():
        pytest.skip(f"{path.name} not available")
    return str(path)


def run(scores, pdb, tmp_path, **kw):
    params = dict(PWES_PARAMS)
    params.update(kw)
    return pwes_clustering(
        df_scores=scores.copy(), pdb_file=pdb,
        out_prefix="test", out_dir=str(tmp_path), **params)


# --------------------------------------------------------------------------- #
# the clustering runs and returns what it promises
# --------------------------------------------------------------------------- #


def test_clustering_returns_six_results(scores, pdb, tmp_path):
    result = run(scores, pdb, tmp_path)
    assert len(result) == 6

    df_gauss, df_sorted, df_unsorted, df_clus, aas_dict, link = result
    assert isinstance(df_gauss, pd.DataFrame) and not df_gauss.empty
    assert isinstance(df_sorted, pd.DataFrame) and not df_sorted.empty
    assert isinstance(df_unsorted, pd.DataFrame) and not df_unsorted.empty
    assert isinstance(df_clus, pd.DataFrame) and not df_clus.empty
    assert isinstance(aas_dict, dict) and len(aas_dict) > 0


def test_pwes_matrix_is_square_and_symmetric(scores, pdb, tmp_path):
    """
    PWES is a pairwise score, so the matrix must be square and symmetric.

    An asymmetry would mean the guide-by-guide product and the distance
    weighting had been combined in the wrong order.
    """
    _, _, df_unsorted, _, _, _ = run(scores, pdb, tmp_path)
    values = df_unsorted.to_numpy(dtype=float)
    assert values.shape[0] == values.shape[1], "PWES matrix is not square"
    np.testing.assert_allclose(values, values.T, rtol=1e-9, atol=1e-9)


def test_gaussian_weights_are_bounded(scores, pdb, tmp_path):
    """
    The distance weight is a Gaussian, so it lies in (0, 1].

    Self-distance is zero, which weights 1; everything else is smaller.
    """
    df_gauss, _, _, _, _, _ = run(scores, pdb, tmp_path)
    values = df_gauss.to_numpy(dtype=float)
    assert np.all(values > 0), "a Gaussian weight cannot be zero or negative"
    assert np.all(values <= 1 + 1e-9), "a Gaussian weight cannot exceed 1"
    np.testing.assert_allclose(np.diag(values), 1.0, atol=1e-9)


def test_every_guide_lands_in_exactly_one_cluster(scores, pdb, tmp_path):
    _, _, _, df_clus, aas_dict, _ = run(scores, pdb, tmp_path)
    assert "cl_new" in df_clus.columns
    counts = df_clus["cl_new"].value_counts()
    assert counts.sum() == len(df_clus)
    assert set(aas_dict) == set(df_clus["cl_new"].unique())


def test_dendrogram_cut_controls_cluster_count(scores, pdb, tmp_path):
    """A tighter cut must not produce fewer clusters than a looser one."""
    _, _, _, loose, _, _ = run(scores, pdb, tmp_path / "loose", dend_t=30)
    _, _, _, tight, _, _ = run(scores, pdb, tmp_path / "tight", dend_t=5)
    assert tight["cl_new"].nunique() >= loose["cl_new"].nunique()


def test_pos_only_keeps_non_negative_scores(scores, pdb, tmp_path):
    """pos_only=True should drop the negative half of the PWES range."""
    _, _, both, _, _, _ = run(scores, pdb, tmp_path / "both", pos_only=False)
    _, _, pos, _, _, _ = run(scores, pdb, tmp_path / "pos", pos_only=True)
    assert both.to_numpy(dtype=float).min() < 0
    assert pos.to_numpy(dtype=float).min() >= -1e-9


# --------------------------------------------------------------------------- #
# outputs land where they are asked to
# --------------------------------------------------------------------------- #


EXPECTED_FILES = [
    "test_aas_dict.pickle", "test_df_gauss.csv", "test_df_pairwise.csv",
    "test_df_pwes_sorted.csv", "test_df_pwes_unsorted.csv", "test_df_clus.csv",
]


def test_outputs_are_written_inside_out_dir(scores, pdb, tmp_path):
    """
    out_dir was concatenated with out_prefix rather than joined.

    Without a trailing separator that created the directory and then wrote
    every file next to it, so a caller passing out_dir='results' ended up with
    an empty results/ and a pile of resultstest_*.csv beside it.
    """
    target = tmp_path / "results"
    run(scores, pdb, target)

    for name in EXPECTED_FILES:
        assert (target / name).exists(), f"{name} was not written into {target}"
    assert not list(tmp_path.glob("resultstest_*")), (
        "files were written beside out_dir instead of inside it"
    )


def test_trailing_separator_in_out_dir_still_works(scores, pdb, tmp_path):
    """The documented usage passed out_dir with a trailing slash."""
    target = tmp_path / "withslash"
    pwes_clustering(df_scores=scores.copy(), pdb_file=pdb,
                    out_prefix="test", out_dir=str(target) + "/", **PWES_PARAMS)
    for name in EXPECTED_FILES:
        assert (target / name).exists(), f"{name} missing"


def test_pickled_clusters_round_trip(scores, pdb, tmp_path):
    run(scores, pdb, tmp_path)
    with open(tmp_path / "test_aas_dict.pickle", "rb") as handle:
        loaded = pickle.load(handle)
    assert isinstance(loaded, dict) and loaded


def test_default_out_dir_does_not_raise(scores, pdb, tmp_path, monkeypatch):
    """
    out_dir defaults to None, which went straight into os.path.exists.

    Every call that did not pass out_dir raised TypeError.
    """
    monkeypatch.chdir(tmp_path)
    pwes_clustering(df_scores=scores.copy(), pdb_file=pdb,
                    out_prefix="test", **PWES_PARAMS)
    assert (tmp_path / "test_df_clus.csv").exists()


def test_caller_frame_is_not_mutated(scores, pdb, tmp_path):
    """pwes_clustering adds an sgRNA_ID column and filters rows."""
    before = scores.copy()
    run(scores, pdb, tmp_path)
    pd.testing.assert_frame_equal(scores, before)


# --------------------------------------------------------------------------- #
# the plots built from the clustering
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def clustered(scores, tmp_path_factory, request):
    path = request.config.rootpath / AF_PDB
    if not path.exists():
        pytest.skip("AF_Q15910.pdb not available")
    out = tmp_path_factory.mktemp("pwes")
    result = pwes_clustering(
        df_scores=scores.copy(), pdb_file=str(path),
        out_prefix="test", out_dir=str(out), **PWES_PARAMS)
    return result, out


def test_pwes_plots_render(clustered):
    from be_scan.figure_plot import (
        pwes_boxplots, pwes_heatmap, pwes_histogram, pwes_scatterplots,
    )
    from be_scan.figure_plot.figure_classes import OutputOpts

    (df_gauss, df_sorted, _, df_clus, _, _), out = clustered
    colors = ["#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#bebada"]
    opts = OutputOpts(save=False, show=False, path=str(out / "fig"))

    pwes_histogram(df_clus, colors, output=opts)
    pwes_boxplots(df_clus, "xpos", colors, output=opts)
    pwes_scatterplots(df_clus, "xpos", "DMSO-PRC2i", colors, output=opts)
    pwes_heatmap(df_sorted, output=opts)


def test_interactive_pwes_plots_render(clustered, tmp_path):
    """The interactive twins must accept the same clustering output."""
    from be_scan.figure_plot import (
        interactive_pwes_boxplots, interactive_pwes_heatmap,
        interactive_pwes_histogram, interactive_pwes_scatterplots,
    )
    from be_scan.figure_plot.figure_classes import OutputOpts

    (df_gauss, df_sorted, _, df_clus, _, _), _ = clustered
    colors = ["#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#bebada"]
    opts = OutputOpts(save=True, show=False, path=str(tmp_path / "fig"))

    interactive_pwes_histogram(df_clus, colors, output=opts)
    interactive_pwes_boxplots(df_clus, "xpos", colors, output=opts)
    interactive_pwes_scatterplots(df_clus, "xpos", "DMSO-PRC2i", colors, output=opts)
    interactive_pwes_heatmap(df_sorted, output=opts)

    written = sorted(p.name for p in tmp_path.glob("*.html"))
    assert len(written) == 4, written
