# The trainer's weighted fits: scikit-learn's bin thresholds at once, the same models

**Date:** 2026-10-12.

**Code:** `audit-fixes` after `c0177ea`. The fix is in the commit that adds this report: `scripts/hgb_binning.py`,
`machine_learning_cmdline.py` (`make_model` installs it; a line in the report header says what it did), and
`HgbBinningTest` in `scripts/test_model_pmml.py`.

**Found in:** [2026-10-12-test-suite-speed](2026-10-12-test-suite-speed/README.md), section 1.

**Machine:** WSL on the shared laptop, CPUs 0-3, on battery. Python 3.14.7 with scikit-learn 1.9.1 and numpy 2.5.3
(`protal-db-build`). OpenMP on one thread.

## The cost

The trainer weighs the scenarios' rows (`--scenario-weight`, 0.25 by default in GTDB builds) and the twins
(`--twin-weight`), so most of its gradient-boosting fits are weighted. scikit-learn 1.8 and later bins a weighted fit
in `sklearn/ensemble/_hist_gradient_boosting/binning.py::_find_binning_thresholds`. For a feature with more than 255
distinct values, it calls `sklearn.utils.stats._weighted_percentile` once per threshold, 254 times. Each call sorts the
column again and goes through the array-API layer. A cProfile of one `--features auto` training on
`TrainerScenariosTest`'s 720 rows: 395 of 402 s in 1.86 million such calls.

One fit, 90 features, 30 rounds, 63 leaves, one thread (the session's `weighted_binning.py`, on a loaded machine):

| rows | unweighted | weighted, scikit-learn's binning |
|---|---|---|
| 2,000 | 1.2 s | 25 s |
| 20,000 | 5.8 s | 53 s |
| 200,000 | 35 s | 268 s |

Above 200,000 rows scikit-learn bins a subsample drawn by the weights, unweighted, so larger fits do not pay this. The
r226 v23 tables have 40,000 to 160,000 rows (`training_data*.tsv` of the test collection), so their fits do: those of
the folds, of every set `--features auto` tries, and the final model.

## The fix

`hgb_binning.py` holds scikit-learn 1.9.1's `_find_binning_thresholds`, but takes the weighted percentiles all at once
in numpy (`_weighted_percentiles`). What the binning passes, a column without NaNs and without weights of 0, gets
scikit-learn's operations in its order:
- the unstable sort of the `(n, 1)` array;
- the cumulative sum of the sorted weights;
- `rank / 100 * total`;
- `searchsorted` on the left, clipped;
- the test against float64's eps, and the average with the next value.

So the thresholds are the same values.

`install()` replaces scikit-learn's function only if, on the columns of `check_columns`, the two agree bit for bit.
Those columns cover:
- many values, ties, and few values;
- NaNs;
- weights of 0, of 1, whole numbers, and 0.25;
- a constant column;
- unweighted columns;
- 255 and 16 bins.

If they differ or anything fails, scikit-learn keeps its own binning. A scikit-learn whose binning takes no weights
(before 1.8, as CI's 1.4.1) is left as it is, and the report says so. The check takes under a second, once per process.
`make_model` calls `install()` before every boosted model, so the worker processes that fit the folds (`--fold-jobs`)
install it too. It imports the module there, and without it they keep scikit-learn's binning: slow, but the same.

## Checks

`HgbBinningTest`:
- `install()` replaces the function, once per process;
- on 12 columns × 2 bin counts (3,000 and 5,000 rows), the thresholds equal the original's, dtype and bits;
- a weighted fit of 2,000 rows (with NaNs and a few-valued feature) gets the same bin thresholds and the same
  probabilities, bit for bit.

Also run: `FoldJobsTest` (the folds in worker processes, weighted) and
`TrainerScenariosTest.test_boosting_by_default_and_the_forest_on_request`. The latter is a weighted training whose
PMML scores as scikit-learn's: `pmml_vs_sklearn_max_diff` 0.0. All pass, 14 s.

One weighted fit, 90 features, the same settings, at load average 0.7:

| rows | unweighted | weighted, scikit-learn's binning | weighted, hgb_binning |
|---|---|---|---|
| 2,000 | 0.46 s | 4.82 s | 0.45 s, the same probabilities |
| 20,000 | 1.05 s | (53 s in the loaded run above) | 1.02 s |

A weighted fit now takes what an unweighted one does.

## Not done (the laptop was on battery)

- A trainer run on a real table before and after. Each fit's binning shrinks as above; how much of a whole training
  that is depends on how much of it is fitting.
- A trainer run as a script with `--fold-jobs` above 1, confirming that its workers import `hgb_binning` (loky passes
  them the parent's `sys.path`). `FoldJobsTest` covers the trainer imported as a module. Either way, a worker that
  cannot import it fits correctly, only slower.
- scikit-learn 1.8 and 1.10 are not installed here. `install()`'s check is what stands for them.

The whole module's tests and the GTDB build end to end were not rerun. The models are bit-identical by construction
and by the checks above.

**The website** is not affected.
