"""
Replicate-based significance in log_transform.

For every condition with more than one replicate, log_transform now emits a
two-sided one-sample t-test of each sgRNA's replicate values against zero,
plus its -log10. Conditions with a single replicate get nothing, because a
p-value across one measurement is not defined.

The tests below check the statistics against scipy directly, confirm the
single-replicate case stays clean, and -- most importantly -- confirm that no
pre-existing column changed.
"""

import numpy as np
import pandas as pd
import pytest
from scipy import stats

from be_scan.analysis.log_transform import (
    calc_moderated_pvals,
    calc_pvals,
    fdr_correct,
    fit_variance_prior,
    log_transform,
    neg_log10,
)


# --------------------------------------------------------------------------- #
# the helpers in isolation
# --------------------------------------------------------------------------- #


def test_calc_pvals_matches_scipy():
    rng = np.random.default_rng(0)
    data = pd.DataFrame(rng.normal(size=(50, 3)), columns=list("abc"))

    got = calc_pvals(data)
    _, expected = stats.ttest_1samp(data.to_numpy(dtype=float), 0.0, axis=1)

    np.testing.assert_allclose(got, expected, rtol=1e-12, atol=0)


def test_calc_pvals_is_small_for_a_consistent_shift():
    """Three replicates all near -3 should be clearly distinguishable from 0."""
    data = pd.DataFrame([[-3.0, -3.1, -2.9]], columns=list("abc"))
    assert calc_pvals(data)[0] < 0.01


def test_calc_pvals_is_large_for_noise_around_zero():
    data = pd.DataFrame([[0.1, -0.15, 0.05]], columns=list("abc"))
    assert calc_pvals(data)[0] > 0.3


def test_calc_pvals_is_nan_without_variance():
    """
    Identical replicates carry no information about significance.

    scipy reports p=0 here, because the t-statistic is infinite, which would
    become -log10(p) ~= 308: the most significant value on the plot. That
    matters in practice, since an sgRNA with zero counts in every replicate
    has exactly identical log-fold-changes. NaN is the honest answer.
    """
    for row in ([1.0, 1.0, 1.0], [-2.5, -2.5, -2.5], [0.0, 0.0, 0.0]):
        data = pd.DataFrame([row], columns=list("abc"))
        assert np.isnan(calc_pvals(data)[0]), f"{row} should give NaN"


def test_dropped_out_guides_do_not_score_as_significant():
    """
    A guide absent from every replicate must not top the volcano plot.

    Zero counts everywhere give identical LFCs, which is exactly the
    zero-variance case above. Verified through the full function rather than
    the helper, since this is how it arises in practice.
    """
    counts, sheet = make_inputs({"treated": 3}, n_guides=10)
    counts.loc[0, ["treated-1", "treated-2", "treated-3"]] = 0
    counts.loc[0, "t0"] = 0

    out = log_transform(sample_sheet=sheet, library_counts=counts,
                        controls=["t0"], save=False, return_df=True)

    dropout = out.loc[0, "treated_LFCminusControl_neglogp"]
    assert np.isnan(dropout), (
        f"a guide with no counts anywhere scored -log10(p)={dropout}"
    )


def test_calc_pvals_handles_partially_missing_replicates():
    """One NaN among three replicates still leaves two usable values."""
    data = pd.DataFrame([[2.0, np.nan, 2.5]], columns=list("abc"))
    assert 0 <= calc_pvals(data)[0] <= 1

    # only one usable replicate is not enough
    data = pd.DataFrame([[2.0, np.nan, np.nan]], columns=list("abc"))
    assert np.isnan(calc_pvals(data)[0])


def test_calc_pvals_is_symmetric_in_sign():
    """A two-sided test cannot care about the direction of the effect."""
    up = pd.DataFrame([[2.0, 2.5, 1.5]], columns=list("abc"))
    down = pd.DataFrame([[-2.0, -2.5, -1.5]], columns=list("abc"))
    np.testing.assert_allclose(calc_pvals(up), calc_pvals(down))


@pytest.mark.parametrize("p,expected", [
    (1.0, 0.0),
    (0.1, 1.0),
    (0.01, 2.0),
    (0.05, 1.3010299956639813),
])
def test_neg_log10(p, expected):
    np.testing.assert_allclose(neg_log10(np.array([p]))[0], expected)


def test_neg_log10_handles_zero_and_nan():
    out = neg_log10(np.array([0.0, np.nan, 0.5]))
    assert np.isfinite(out[0]) and out[0] > 300, "p=0 must not become inf"
    assert np.isnan(out[1])
    np.testing.assert_allclose(out[2], -np.log10(0.5))


# --------------------------------------------------------------------------- #
# fixtures for the full function
# --------------------------------------------------------------------------- #


def make_inputs(rep_counts, n_guides=40, seed=0):
    """
    Build a counts table and sample sheet with the given replicate structure.

    rep_counts maps condition name -> number of replicates. A 't0' control
    with a single sample is always included.
    """
    rng = np.random.default_rng(seed)
    counts = {"sgRNA_seq": [f"guide_{i}" for i in range(n_guides)]}
    rows = [{"condition": "t0", "agg_conditions": "t0", "counts_file": "t0.csv"}]
    counts["t0"] = rng.integers(200, 400, n_guides)

    for cond, n in rep_counts.items():
        for r in range(1, n + 1):
            name = f"{cond}-{r}"
            counts[name] = rng.integers(200, 400, n_guides)
            rows.append({"condition": name, "agg_conditions": cond,
                         "counts_file": f"{name}.csv"})

    return pd.DataFrame(counts), pd.DataFrame(rows)


def run(rep_counts, **kw):
    counts, sheet = make_inputs(rep_counts, **kw)
    return log_transform(sample_sheet=sheet, library_counts=counts,
                         controls=["t0"], save=False, return_df=True)


# --------------------------------------------------------------------------- #
# the columns appear, and only where they should
# --------------------------------------------------------------------------- #


def test_pvalue_columns_added_for_multi_replicate_conditions():
    out = run({"treated": 3})
    assert "treated_LFCminusControl_pval" in out.columns
    assert "treated_LFCminusControl_neglogp" in out.columns


def test_no_pvalue_columns_for_single_replicate_conditions():
    """A p-value across one replicate is not defined, so none is produced."""
    out = run({"solo": 1})
    assert "solo_LFCminusControl_pval" not in out.columns
    assert "solo_LFCminusControl_neglogp" not in out.columns
    # the pre-existing columns are still there
    assert "solo_LFCminusControl_avg" in out.columns
    assert "solo_LFCminusControl_stdev" in out.columns


def test_mixed_replicate_structure():
    """Multi-rep conditions get p-values; single-rep ones alongside do not."""
    out = run({"treated": 3, "solo": 1, "pair": 2})
    for cond in ("treated", "pair"):
        assert f"{cond}_LFCminusControl_pval" in out.columns
    assert "solo_LFCminusControl_pval" not in out.columns


def test_pvalues_are_in_range():
    out = run({"treated": 4})
    p = out["treated_LFCminusControl_pval"].dropna()
    assert len(p) > 0
    assert ((p >= 0) & (p <= 1)).all(), "p-values outside [0,1]"


def test_neglogp_is_the_negative_log10_of_pval():
    out = run({"treated": 3})
    p = out["treated_LFCminusControl_pval"]
    q = out["treated_LFCminusControl_neglogp"]
    mask = p.notna() & (p > 0)
    np.testing.assert_allclose(q[mask], -np.log10(p[mask]), rtol=1e-12)


def test_pvalue_uses_the_same_replicates_as_the_average():
    """The test must run on the same columns that feed _avg and _stdev."""
    counts, sheet = make_inputs({"treated": 3})
    out = log_transform(sample_sheet=sheet, library_counts=counts,
                        controls=["t0"], save=False, return_df=True)

    reps = [f"treated-{r}_LFCminusControl" for r in (1, 2, 3)]
    expected = calc_pvals(out[reps])
    np.testing.assert_allclose(
        out["treated_LFCminusControl_pval"].to_numpy(), expected, rtol=1e-12)


def test_works_without_controls():
    """With controls=[] the suffix is _LFC, and p-values follow it."""
    counts, sheet = make_inputs({"treated": 3})
    out = log_transform(sample_sheet=sheet, library_counts=counts,
                        controls=[], save=False, return_df=True)
    assert "treated_LFC_pval" in out.columns
    assert "treated_LFC_neglogp" in out.columns


# --------------------------------------------------------------------------- #
# nothing that existed before has changed
# --------------------------------------------------------------------------- #


PREEXISTING_SUFFIXES = ["_LFC", "_LFCminusControl"]


def test_existing_columns_are_untouched():
    """
    Every column v2.0.0 produced must still be present, with the same values.

    The p-value work must be purely additive: downstream code reads these
    frames by column name.
    """
    counts, sheet = make_inputs({"treated": 3, "other": 2})
    out = log_transform(sample_sheet=sheet, library_counts=counts,
                        controls=["t0"], save=False, return_df=True)

    expected = [
        "sgRNA_seq", "t0", "treated-1", "treated-2", "treated-3",
        "other-1", "other-2", "control_avg", "control_avg_LFC", "t0_LFC",
        "treated-1_LFC", "treated-1_LFCminusControl",
        "treated_LFCminusControl_avg", "treated_LFCminusControl_stdev",
        "other_LFCminusControl_avg", "other_LFCminusControl_stdev",
    ]
    for col in expected:
        assert col in out.columns, f"{col} disappeared"


def test_new_columns_are_appended_not_interleaved():
    """
    Existing columns keep their relative order.

    Code doing df.iloc[:, :n] or reading positionally still works.
    """
    counts, sheet = make_inputs({"treated": 3})
    out = log_transform(sample_sheet=sheet, library_counts=counts,
                        controls=["t0"], save=False, return_df=True)

    cols = list(out.columns)
    avg = cols.index("treated_LFCminusControl_avg")
    stdev = cols.index("treated_LFCminusControl_stdev")
    pval = cols.index("treated_LFCminusControl_pval")
    neglogp = cols.index("treated_LFCminusControl_neglogp")

    assert avg < stdev < pval < neglogp, (
        "p-value columns must follow the existing avg/stdev columns"
    )


def test_averages_are_unaffected_by_the_new_columns():
    """The _avg and _stdev values must be exactly what they always were."""
    counts, sheet = make_inputs({"treated": 3})
    out = log_transform(sample_sheet=sheet, library_counts=counts,
                        controls=["t0"], save=False, return_df=True)

    reps = [f"treated-{r}_LFCminusControl" for r in (1, 2, 3)]
    np.testing.assert_allclose(
        out["treated_LFCminusControl_avg"], out[reps].mean(axis=1))
    np.testing.assert_allclose(
        out["treated_LFCminusControl_stdev"], out[reps].std(axis=1))


# --------------------------------------------------------------------------- #
# variance shrinkage
# --------------------------------------------------------------------------- #


def test_variance_prior_recovers_a_known_distribution():
    """
    The fitted prior should land near the truth for simulated data.

    Variances are drawn as sigma^2 * chi2(df)/df around a known sigma^2, which
    is exactly the model fit_variance_prior assumes.
    """
    rng = np.random.default_rng(7)
    true_var, df = 0.25, 2
    variances = true_var * rng.chisquare(df, size=4000) / df

    prior_var, prior_df = fit_variance_prior(variances, df)

    assert prior_var == pytest.approx(true_var, rel=0.25), prior_var
    assert prior_df > 1, "a genuinely dispersed set should give finite prior df"


def test_variance_prior_handles_degenerate_input():
    assert np.isnan(fit_variance_prior(np.array([]), 2)[0])
    assert np.isnan(fit_variance_prior(np.array([0.0, 0.0]), 2)[0])


def test_moderated_test_is_defined_without_variance():
    """
    Shrinkage rescues the case the ordinary test cannot handle.

    Identical replicates give the ordinary test NaN. The moderated test
    borrows the prior variance and stays finite.
    """
    rng = np.random.default_rng(3)
    data = pd.DataFrame(rng.normal(scale=0.4, size=(300, 3)))
    data.iloc[0] = [1.0, 1.0, 1.0]   # zero variance

    assert np.isnan(calc_pvals(data)[0])
    assert np.isfinite(calc_moderated_pvals(data)[0])


def test_moderated_test_penalises_implausibly_small_variance():
    """
    A guide with a tiny effect and a lucky-tight variance must not win.

    This is the failure mode that puts non-targeting controls at the top of a
    3-replicate screen: the ordinary test rewards agreement, the moderated
    test asks whether the effect is large relative to typical spread.
    """
    rng = np.random.default_rng(11)
    data = pd.DataFrame(rng.normal(scale=0.5, size=(500, 3)))
    data.iloc[0] = [0.02, 0.021, 0.0205]   # negligible effect, tiny spread
    data.iloc[1] = [1.6, 2.0, 1.2]         # large effect, ordinary spread

    plain = calc_pvals(data)
    moderated = calc_moderated_pvals(data)

    assert plain[0] < plain[1], "ordinary test should favour the tight guide"
    assert moderated[1] < moderated[0], (
        "moderated test should favour the guide with the real effect"
    )


def test_moderated_pvalues_are_in_range():
    rng = np.random.default_rng(5)
    data = pd.DataFrame(rng.normal(size=(200, 3)))
    p = calc_moderated_pvals(data)
    assert np.all((p[np.isfinite(p)] >= 0) & (p[np.isfinite(p)] <= 1))


def test_moderated_test_needs_two_replicates():
    data = pd.DataFrame([[1.0, np.nan, np.nan]] * 50)
    assert np.isnan(calc_moderated_pvals(data)[0])


# --------------------------------------------------------------------------- #
# FDR correction
# --------------------------------------------------------------------------- #


def test_fdr_matches_statsmodels():
    import statsmodels.stats.multitest as smm
    rng = np.random.default_rng(1)
    p = rng.uniform(size=200)
    _, expected, _, _ = smm.multipletests(p, method="fdr_bh")
    np.testing.assert_allclose(fdr_correct(p), expected)


def test_fdr_is_never_smaller_than_the_raw_pvalue():
    rng = np.random.default_rng(2)
    p = rng.uniform(size=500)
    q = fdr_correct(p)
    assert np.all(q >= p - 1e-12), "BH q-values cannot be below their p-values"


def test_fdr_preserves_nan_and_excludes_it_from_correction():
    p = np.array([0.001, np.nan, 0.5, 0.9])
    q = fdr_correct(p)
    assert np.isnan(q[1])
    assert np.isfinite(q[0]) and np.isfinite(q[2]) and np.isfinite(q[3])


def test_fdr_of_all_nan_is_all_nan():
    assert np.all(np.isnan(fdr_correct(np.array([np.nan, np.nan]))))


# --------------------------------------------------------------------------- #
# the moderated and FDR columns in the full function
# --------------------------------------------------------------------------- #


def test_all_six_significance_columns_are_added():
    out = run({"treated": 3})
    for suffix in ("_pval", "_neglogp", "_pval_mod", "_neglogp_mod",
                   "_qval", "_qval_mod"):
        assert f"treated_LFCminusControl{suffix}" in out.columns, suffix


def test_no_significance_columns_for_single_replicates():
    out = run({"solo": 1})
    for suffix in ("_pval", "_neglogp", "_pval_mod", "_neglogp_mod",
                   "_qval", "_qval_mod"):
        assert f"solo_LFCminusControl{suffix}" not in out.columns, suffix


def test_qval_is_the_fdr_of_its_own_pval():
    """Each q-value column must correct the p-value column it is named after."""
    out = run({"treated": 3})
    for raw, corrected in (("_pval", "_qval"), ("_pval_mod", "_qval_mod")):
        p = out[f"treated_LFCminusControl{raw}"].to_numpy()
        q = out[f"treated_LFCminusControl{corrected}"].to_numpy()
        np.testing.assert_allclose(q, fdr_correct(p), equal_nan=True)


def test_saved_csv_round_trips(tmp_path):
    counts, sheet = make_inputs({"treated": 3})
    log_transform(sample_sheet=sheet, library_counts=counts, controls=["t0"],
                  out_dir=str(tmp_path), out_file="lfc.csv",
                  save=True, return_df=False)

    written = pd.read_csv(tmp_path / "lfc.csv")
    assert "treated_LFCminusControl_pval" in written.columns
    assert "treated_LFCminusControl_neglogp" in written.columns
