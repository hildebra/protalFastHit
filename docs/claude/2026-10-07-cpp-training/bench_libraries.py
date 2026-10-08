#!/usr/bin/env python3
"""Fit time and quality of protal's presence model in scikit-learn (the trainer's) and in the C++ tree libraries
LightGBM and XGBoost, on the same training table, features, row weights and settings.

The trainer's boosted model (machine_learning_cmdline.py, defaults since 15006c2): 250 rounds at a learning rate of
0.1, at most 63 leaves, at least 20 rows a leaf, L2 1, 255 bins, classes balanced, the scenarios' rows at weight 0.25.
LightGBM takes the same settings by name. XGBoost (hist, leaf-wise) has no minimum of rows per leaf: its
min_child_weight (a minimum of hessian) stays at its default, 1.

    python3 bench_libraries.py --truth-file train_pe.tsv --threads 1,2,4 --repeats 2 --out fit_times.tsv
    python3 bench_libraries.py --truth-file train_pe.tsv --threads 4 --folds 5 --out quality.tsv   # species held out
    python3 bench_libraries.py --truth-file train_pe.tsv --threads 4 --forest --out forest.tsv     # --model forest

Written: one row per library, threads and repeat (fit and predict seconds), and with --folds the scores with species
held out (GroupKFold on the taxon, as the trainer's "species" rows): F1 at 0.5, AP, AUC, log loss.
"""

import argparse
import os
import sys
import time

import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import average_precision_score, f1_score, log_loss, roc_auc_score
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

sys.path.insert(0, os.environ.get("PROTAL_SCRIPTS", os.path.join(os.path.dirname(__file__), "..", "..", "..", "scripts")))
import machine_learning_cmdline as mlc  # noqa: E402
from model_features import feature_columns  # noqa: E402

FEATURES = "normalized+adjacency+distance+depth+divergence+unfiltered+ref"  # the r226 v14 build's set
ROUNDS, RATE, LEAVES, MIN_LEAF, L2, BINS = 250, 0.1, 63, 20, 1.0, 255
FOREST_TREES, FOREST_LEAVES = 64, 256


def balanced(y):
    """scikit-learn's class_weight="balanced" as HistGradientBoostingClassifier applies it: n / (2 * rows of the class)."""
    counts = np.bincount(y, minlength=2)
    return (len(y) / (2.0 * counts))[y]


def gbm(library, threads):
    if library == "sklearn":
        return HistGradientBoostingClassifier(max_iter=ROUNDS, learning_rate=RATE, max_leaf_nodes=LEAVES,
                                              min_samples_leaf=MIN_LEAF, l2_regularization=L2, max_bins=BINS,
                                              class_weight="balanced", early_stopping=False, random_state=1)
    if library == "lightgbm":
        import lightgbm
        return lightgbm.LGBMClassifier(n_estimators=ROUNDS, learning_rate=RATE, num_leaves=LEAVES, max_depth=-1,
                                       min_child_samples=MIN_LEAF, min_child_weight=1e-3, reg_lambda=L2, max_bin=BINS,
                                       n_jobs=threads, random_state=1, verbose=-1)
    if library == "xgboost":
        import xgboost
        return xgboost.XGBClassifier(n_estimators=ROUNDS, learning_rate=RATE, tree_method="hist", grow_policy="lossguide",
                                     max_leaves=LEAVES, max_depth=0, reg_lambda=L2, max_bin=BINS, n_jobs=threads,
                                     random_state=1)
    raise ValueError(library)


def forest(library, threads, n_features):
    """--model forest: 64 trees of at most 256 leaves, sqrt of the features at each split, bootstrap rows."""
    share = np.sqrt(n_features) / n_features
    if library == "sklearn":
        return RandomForestClassifier(n_estimators=FOREST_TREES, max_leaf_nodes=FOREST_LEAVES, min_samples_leaf=1,
                                      max_features="sqrt", class_weight="balanced", n_jobs=threads, random_state=1)
    if library == "lightgbm":
        import lightgbm
        return lightgbm.LGBMClassifier(boosting_type="rf", n_estimators=FOREST_TREES, num_leaves=FOREST_LEAVES,
                                       min_child_samples=1, min_child_weight=1e-3, bagging_fraction=0.632,
                                       bagging_freq=1, feature_fraction_bynode=share, max_bin=BINS, n_jobs=threads,
                                       random_state=1, verbose=-1)
    if library == "xgboost":
        import xgboost
        return xgboost.XGBRFClassifier(n_estimators=FOREST_TREES, tree_method="hist", grow_policy="lossguide",
                                       max_leaves=FOREST_LEAVES, max_depth=0, subsample=0.632, colsample_bynode=share,
                                       reg_lambda=0, min_child_weight=1e-3, max_bin=BINS, n_jobs=threads, random_state=1)
    raise ValueError(library)


def make(library, threads, n_features, is_forest):
    return forest(library, threads, n_features) if is_forest else gbm(library, threads)


def fit(model, library, X, y, w):
    """The fit with the trainer's weights: the scenarios' row weights, and for the libraries other than scikit-learn
    (whose estimator balances the classes itself) the balanced class weights multiplied in."""
    if library == "sklearn":
        return model.fit(X, y, sample_weight=w)
    return model.fit(X, y, sample_weight=w * balanced(y))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--truth-file", required=True)
    p.add_argument("--libraries", default="sklearn,lightgbm,xgboost")
    p.add_argument("--threads", default="4", help="comma-separated thread counts")
    p.add_argument("--repeats", type=int, default=1)
    p.add_argument("--folds", type=int, default=0, help="also score with species held out in this many folds")
    p.add_argument("--forest", action="store_true", help="the random forest (--model forest) instead of boosting")
    p.add_argument("--rows", type=int, default=0,
                   help="fit on this many rows drawn at random (e.g. 4000: the fixed cost of a fit's 250 trees, whose "
                        "per-row work is then small)")
    p.add_argument("--out", required=True)
    opts = p.parse_args()

    t0 = time.time()
    df = mlc.load_table(opts.truth_file)
    cols = [c for c in feature_columns(df.columns, FEATURES) if c != "domain"]
    X = np.ascontiguousarray(df[cols].to_numpy(dtype=np.float64))
    y = df["truth"].to_numpy().astype(int)
    w = np.where(mlc.scenario_of(df) != "", 0.25, 1.0)
    groups = df["taxon"].astype(str).to_numpy()
    if opts.rows:
        keep = np.sort(np.random.default_rng(1).choice(len(y), size=min(opts.rows, len(y)), replace=False))
        X, y, w, groups = np.ascontiguousarray(X[keep]), y[keep], w[keep], groups[keep]
    print(f"{len(y)} rows, {len(cols)} features, loaded in {time.time() - t0:.1f} s", flush=True)

    libraries = opts.libraries.split(",")
    rows = []
    for threads in [int(t) for t in opts.threads.split(",")]:
        for repeat in range(opts.repeats):
            # alternate the libraries' order between repeats, so that a slow spell of the machine hits each
            order = libraries if repeat % 2 == 0 else libraries[::-1]
            for library in order:
                with threadpool_limits(limits=threads, user_api="openmp"):
                    model = make(library, threads, len(cols), opts.forest)
                    t = time.time()
                    fit(model, library, X, y, w)
                    fit_s = time.time() - t
                    t = time.time()
                    prob = model.predict_proba(X)[:, 1]
                    predict_s = time.time() - t
                row = dict(library=library, threads=threads, repeat=repeat, fit_s=round(fit_s, 2),
                           predict_s=round(predict_s, 2), train_log_loss=round(log_loss(y, prob), 5))
                print(row, flush=True)
                rows.append(row)

    if opts.folds:
        threads = int(opts.threads.split(",")[-1])
        splits = list(GroupKFold(opts.folds).split(X, y, groups))
        for library in libraries:
            prob = np.full(len(y), np.nan)
            t = time.time()
            with threadpool_limits(limits=threads, user_api="openmp"):
                for train, test in splits:
                    model = fit(make(library, threads, len(cols), opts.forest), library, X[train], y[train], w[train])
                    prob[test] = model.predict_proba(X[test])[:, 1]
            row = dict(library=library, threads=threads, repeat="species_folds", fit_s=round(time.time() - t, 2),
                       predict_s="", train_log_loss="",
                       F1=round(f1_score(y, prob >= 0.5), 4), AP=round(average_precision_score(y, prob), 4),
                       AUC=round(roc_auc_score(y, prob), 4), log_loss=round(log_loss(y, prob), 4))
            print(row, flush=True)
            rows.append(row)

    keys = ["library", "threads", "repeat", "fit_s", "predict_s", "train_log_loss", "F1", "AP", "AUC", "log_loss"]
    with open(opts.out, "w") as out:
        out.write("\t".join(keys) + "\n")
        for row in rows:
            out.write("\t".join(str(row.get(k, "")) for k in keys) + "\n")


if __name__ == "__main__":
    main()
