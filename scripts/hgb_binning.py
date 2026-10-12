#!/usr/bin/env python3
"""hgb_binning.py - the bin thresholds of a weighted HistGradientBoostingClassifier fit, as scikit-learn computes them,
in one pass per feature instead of one per threshold: the same thresholds, so the same models.

scikit-learn 1.8 and later (sklearn/ensemble/_hist_gradient_boosting/binning.py, _find_binning_thresholds) bins a feature
of a weighted fit that has more than max_bins (255) distinct values by calling sklearn.utils.stats._weighted_percentile
once per threshold, 254 times, each call sorting the column again. One fit of 90 features and 30 rounds on one thread
took 25 s instead of 1.2 s on 2,000 rows, and 268 s instead of 35 s on 200,000 (above 200,000 rows scikit-learn bins a
weighted subsample, unweighted). The trainer weighs the scenarios' rows (--scenario-weight) and the twins
(--twin-weight), so the fits of every GTDB build's training paid this (docs/claude/2026-10-12-trainer-weighted-binning).
Here the 254 percentiles are computed at once in numpy (_weighted_percentiles), with scikit-learn's operations in its
order.

install() puts _find_binning_thresholds below, scikit-learn 1.9.1's with that change, in place of scikit-learn's, if
the two give the same thresholds, bit for bit, on columns of every kind (check_columns: many values or few, ties, NaNs,
weights of 0, constant, unweighted); otherwise scikit-learn is left as it is (scikit-learn before 1.8 bins without the
weights). machine_learning_cmdline.make_model calls it, so it is in place in the processes that fit the folds too.

  python3 -m unittest scripts/test_model_pmml.py     # HgbBinningTest
"""

import numpy as np

_installed = None  # install()'s answer, once it has run in this process
_original = None   # scikit-learn's _find_binning_thresholds, once install() has replaced it


def _weighted_percentiles(array, sample_weight, percentile_rank):
    """scikit-learn 1.9.1's sklearn.utils.stats._weighted_percentile(array, sample_weight, percentile_rank, average=True)
    for what the binning passes it, a 1-D array without NaNs and weights without zeros: every rank at once in numpy,
    where scikit-learn takes one rank after another through the array API. The same operations in the same order (its
    sort of the (n, 1) array, the cumulative sum, rank / 100 * total, searchsorted on the left, the eps test and the
    average of the two values), so the same values."""
    array = np.asarray(array, dtype=np.float64)
    sample_weight = np.asarray(sample_weight, dtype=np.float64)
    percentile_rank = np.asarray(percentile_rank, dtype=np.float64)
    last = array.shape[0] - 1
    sorted_idx = np.argsort(array.reshape(-1, 1), axis=0)[:, 0]
    weight_cdf = np.cumsum(sample_weight[sorted_idx])
    adjusted = percentile_rank / 100 * weight_cdf[-1]
    zero = adjusted == 0  # ignore leading weights of 0 at rank 0, as scikit-learn does
    adjusted[zero] = np.nextafter(adjusted[zero], adjusted[zero] + 1)
    indices = np.clip(np.searchsorted(weight_cdf, adjusted), 0, last)
    in_sorted = sorted_idx[indices]
    fraction_above = weight_cdf[indices] - adjusted
    plus_one_in_sorted = sorted_idx[np.clip(indices + 1, 0, last)]
    return np.where(fraction_above > np.finfo(np.float64).eps, array[in_sorted],
                    (array[in_sorted] + array[plus_one_in_sorted]) / 2)


def _find_binning_thresholds(col_data, max_bins, sample_weight=None):
    """scikit-learn 1.9.1's sklearn.ensemble._hist_gradient_boosting.binning._find_binning_thresholds (BSD-3-Clause, the
    scikit-learn developers), with the weighted percentiles computed all at once (_weighted_percentiles)."""
    from sklearn.ensemble._hist_gradient_boosting.common import ALMOST_INF, X_DTYPE
    from numpy.lib.stride_tricks import sliding_window_view

    # ignore missing values when computing bin thresholds
    missing_mask = np.isnan(col_data)
    any_missing = missing_mask.any()
    if any_missing:
        col_data = col_data[~missing_mask]

    # If sample_weight is not None and 0-weighted values exist, we need to
    # remove those before calculating the distinct points.
    if sample_weight is not None:
        if any_missing:
            sample_weight = sample_weight[~missing_mask]
        nnz_sw = sample_weight != 0
        col_data = col_data[nnz_sw]
        sample_weight = sample_weight[nnz_sw]

    # The data will be sorted anyway in np.unique and again in percentile, so we do it
    # here. Sorting also returns a contiguous array.
    sort_idx = np.argsort(col_data)
    col_data = col_data[sort_idx]
    if sample_weight is not None:
        sample_weight = sample_weight[sort_idx]

    distinct_values = np.unique(col_data).astype(X_DTYPE)

    if len(distinct_values) == 1:
        return np.asarray([])

    if len(distinct_values) <= max_bins:
        # Calculate midpoints if distinct values <= max_bins
        bin_thresholds = sliding_window_view(distinct_values, 2).mean(axis=1)
    elif sample_weight is None:
        # We compute bin edges using the output of np.percentile with
        # the "averaged_inverted_cdf" interpolation method that is consistent
        # with the code for the sample_weight != None case.
        percentiles = np.linspace(0, 100, num=max_bins + 1)
        percentiles = percentiles[1:-1]
        bin_thresholds = np.percentile(
            col_data, percentiles, method="averaged_inverted_cdf"
        )
        assert bin_thresholds.shape[0] == max_bins - 1
    else:
        percentiles = np.linspace(0, 100, num=max_bins + 1)
        percentiles = percentiles[1:-1]
        # every percentile at once (scikit-learn: _weighted_percentile once for each), from the column without its NaNs
        # and weights of 0 (removed above)
        bin_thresholds = _weighted_percentiles(col_data, sample_weight, percentiles)
        assert bin_thresholds.shape[0] == max_bins - 1
    # Remove duplicated thresholds if they exist.
    unique_bin_values = np.unique(bin_thresholds)
    if unique_bin_values.shape[0] != bin_thresholds.shape[0]:
        bin_thresholds = unique_bin_values

    # We avoid having +inf thresholds: +inf thresholds are only allowed in
    # a "split on nan" situation.
    np.clip(bin_thresholds, a_min=None, a_max=ALMOST_INF, out=bin_thresholds)
    return bin_thresholds


def check_columns(n=400, seed=7):
    """[(name, column, weights)] on which install() compares the two functions: each branch of the binning (many values
    or 255 and fewer, constant, unweighted), with ties, NaNs, weights of 0 and whole-number weights."""
    rng = np.random.default_rng(seed)
    weights = np.where(rng.random(n) < 0.3, 0.25, 1.0)
    zero = weights.copy()
    zero[rng.random(n) < 0.2] = 0.0
    many = rng.normal(size=n)
    ties = np.round(rng.normal(size=n) * 100) / 7  # ~300 values, many of them twice or more
    nans = many.copy()
    nans[rng.random(n) < 0.1] = np.nan
    return [("many values", many, weights), ("ties", ties, weights), ("NaNs", nans, weights),
            ("weights of 0", many, zero), ("ties and weights of 0", ties, zero),
            ("whole-number weights", many, rng.integers(1, 4, n).astype(float)), ("weights of 1", ties, np.ones(n)),
            ("few values", np.round(many, 1), weights), ("constant", np.ones(n), weights),
            ("unweighted", many, None), ("unweighted, few values", np.round(many, 1), None)]


def install():
    """Put _find_binning_thresholds in place of scikit-learn's, if both give the same thresholds on check_columns, bit
    for bit; once per process -> what it did, in words (for the trainer's report)."""
    global _installed, _original
    if _installed is not None:
        return _installed
    try:
        import inspect
        import sklearn
        from sklearn.ensemble._hist_gradient_boosting import binning
        original = binning._find_binning_thresholds
        if original is _find_binning_thresholds:
            _installed = "weighted bin thresholds in one pass per feature (hgb_binning.py)"
            return _installed
        if "sample_weight" not in inspect.signature(original).parameters:
            _installed = f"scikit-learn {sklearn.__version__} bins without the sample weights: its own binning"
            return _installed
        for name, column, weights in check_columns():
            for max_bins in (255, 16):
                theirs = original(column.copy(), max_bins, sample_weight=None if weights is None else weights.copy())
                ours = _find_binning_thresholds(column.copy(), max_bins,
                                                sample_weight=None if weights is None else weights.copy())
                if theirs.dtype != ours.dtype or not np.array_equal(theirs, ours):
                    raise ValueError(f"other thresholds on the column with {name}, max_bins {max_bins}")
    except Exception as error:  # an older or newer scikit-learn: its own binning, slow or not
        _installed = f"scikit-learn's own binning ({type(error).__name__}: {error})"
        return _installed
    _original, binning._find_binning_thresholds = original, _find_binning_thresholds
    _installed = (f"weighted bin thresholds in one pass per feature (hgb_binning.py; scikit-learn {sklearn.__version__} "
                  "computes them one at a time)")
    return _installed
