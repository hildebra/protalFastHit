#!/usr/bin/env python3
"""What the new feature groups add to the r226 v17 models: refits on v17's tables without them, species held out
(5 folds grouped by taxon) and on the test set (the final model), as the 2026-10-07 v15 report's ablate.py did but
without the sample folds.

Variants (gradient boosting as v17 trained it: 250 rounds at 0.1, 63 leaves, min leaf 20, L2 1, balanced classes;
scenario rows weighted 0.25):
  v17           the build's 81 features
  no_ancestry   without ancestry_sites_per_record, ancestry_agreement, ancestry_congener_share
  no_foreign    without foreign_scanned_share, foreign_copy_share, foreign_genus_copy_share
  no_gaps       without gap_*, untried_candidate_rate
  no_new        without ancestry, foreign, gaps, untried (the 70 of 2026-10-07)
  v15_set       the v15 build's 55 features (without the consistency, shape and neighbourhood groups too)

Writes OUT/<read type>_<variant>.tsv.gz (training rows with p_species, test rows with p) and OUT/results.tsv: per
read type and variant the knob chosen with species held out, F1 with species held out at 0.5 and at the knob, and
on the test set at 0.5, at the knob and at the best test threshold, with FP and FN.

    python3 ablate_v17.py --build local/v17 --out ~/v17/ablate --read-types pe -t 3
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import machine_learning_cmdline as trainer  # noqa: E402
import model_features as mf  # noqa: E402

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
ANCESTRY = ["ancestry_sites_per_record", "ancestry_agreement", "ancestry_congener_share"]
FOREIGN = ["foreign_scanned_share", "foreign_copy_share", "foreign_genus_copy_share"]
GAPS = ["gap_informative_share", "gap_within_min_share", "gap_within_median_share", "gap_position", "untried_candidate_rate"]
V15_SET = "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity"
VARIANTS = ["v17", "no_ancestry", "no_foreign", "no_gaps", "no_new", "v15_set"]
KEEP = ["meta_sample", "meta_read_pairs", "meta_scenario", "taxon", "truth", "identity", "fragments"]


def columns_of(variant, columns):
    full = [c for c in mf.feature_set_columns(mf.DEFAULT_FEATURE_SET) if c in columns]
    extra = [c for c in FOREIGN + GAPS if c in columns and c not in full]
    cols = full + extra
    if variant == "v15_set":
        return [c for c in mf.feature_set_columns(V15_SET) if c in columns]
    drop = {"v17": [], "no_ancestry": ANCESTRY, "no_foreign": FOREIGN, "no_gaps": GAPS,
            "no_new": ANCESTRY + FOREIGN + GAPS}[variant]
    return [c for c in cols if c not in drop]


def model(seed):
    return HistGradientBoostingClassifier(max_iter=250, learning_rate=0.1, max_leaf_nodes=63, min_samples_leaf=20,
                                          l2_regularization=1.0, class_weight="balanced", early_stopping=False,
                                          random_state=seed)


def counts(y, p, t, w=None):
    w = np.ones(len(y)) if w is None else w
    call = p >= t
    tp, fp, fn = (w * (call & (y == 1))).sum(), (w * (call & (y == 0))).sum(), (w * (~call & (y == 1))).sum()
    return tp, fp, fn, (2 * tp / (2 * tp + fp + fn) if tp else 0.0)


def choose_knob(y, p, w):
    grid = np.round(np.arange(0.30, 0.851, 0.01), 2)
    f = np.array([counts(y, p, t, w)[3] for t in grid])
    return float(grid[f.argmax()]) if f.max() >= counts(y, p, 0.5, w)[3] + 0.002 else 0.5


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--read-types", default="pe")
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("-t", "--threads", type=int, default=3)
    opts = ap.parse_args()
    warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
    os.makedirs(opts.out, exist_ok=True)
    results = os.path.join(opts.out, "results.tsv")
    for rt in opts.read_types.split(","):
        train = trainer.load_table(os.path.join(opts.build, "training", TABLES[rt]))
        test = trainer.load_table(os.path.join(opts.build, "test", TABLES[rt]))
        y, yt = train["truth"].to_numpy(), test["truth"].to_numpy()
        sc = trainer.scenario_of(train)
        w = np.where(sc == "", 1.0, 0.25)
        species = list(GroupKFold(opts.folds).split(train, y, train["taxon"].astype(str)))
        for variant in opts.variants.split(","):
            began = time.time()
            cols = columns_of(variant, train.columns)
            X = train[cols].to_numpy(dtype=np.float64)
            Xt = test[cols].to_numpy(dtype=np.float64)
            p_species = np.full(len(y), np.nan)
            with threadpool_limits(limits=opts.threads, user_api="openmp"):
                for tr, te in species:
                    p_species[te] = model(opts.seed).fit(X[tr], y[tr], sample_weight=w[tr]).predict_proba(X[te])[:, 1]
                p_test = model(opts.seed).fit(X, y, sample_weight=w).predict_proba(Xt)[:, 1]
            knob = choose_knob(y, p_species, w)
            grid = np.round(np.arange(0.05, 0.96, 0.01), 2)
            best = max(grid, key=lambda t: counts(yt, p_test, t)[3])
            row = {"read type": rt, "variant": variant, "features": len(cols), "knob": knob,
                   "held-out F1 at 0.5": counts(y, p_species, 0.5)[3], "held-out F1 at knob": counts(y, p_species, knob)[3],
                   "test F1 at 0.5": counts(yt, p_test, 0.5)[3], "test F1 at knob": counts(yt, p_test, knob)[3],
                   "test FP/FN at knob": "%d/%d" % counts(yt, p_test, knob)[1:3],
                   "best test F1": counts(yt, p_test, best)[3], "at": best, "seconds": round(time.time() - began)}
            a = train[KEEP].copy()
            a.insert(0, "part", "training")
            a["p_species"], a["p"] = p_species, np.nan
            b = test[KEEP].copy()
            b.insert(0, "part", "test")
            b["p_species"], b["p"] = np.nan, p_test
            pd.concat([a, b]).to_csv(os.path.join(opts.out, f"{rt}_{variant}.tsv.gz"), sep="\t", index=False,
                                     float_format="%.6g")
            pd.DataFrame([row]).round(4).to_csv(results, sep="\t", index=False, mode="a", header=not os.path.exists(results))
            print(" ".join(f"{k}={v}" for k, v in row.items()), flush=True)


if __name__ == "__main__":
    main()
