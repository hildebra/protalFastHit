#!/usr/bin/env python3
"""What did the sample's complexity features add to the r226 v15 models? Refits on v15's tables, with row predictions
kept for paired comparisons (bootstrap.py).

Variants (all gradient boosting as v15 trained it: 250 rounds at 0.1, 63 leaves, min leaf 20, L2 1, balanced
classes, no early stopping; scenario rows weighted 0.25):

  v15            normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity (the build's 55 features)
  no_complexity  without the complexity group (sample_log_taxa, sample_low_identity, sample_identity): v14's features
  no_depth       v15 without sample_log_fragments
  no_ref         v15 without the ref group (su_rate_ref, lu_rate_ref, lsu_rate_ref: constant within a species, they
                 can name it, and the final model can learn which species tend to be present)

Each variant: 5 folds grouped by taxon (species held out: p_species, and the knob as the trainer chooses its global
one: the highest weighted F1 over 0.30-0.85, kept if it gains 0.002 over 0.5), 5 folds grouped by sample (samples
held out: p_samples), and the final model on all rows, which scores the test set (the design's test samples and the
scenarios' hold-out samples; p).

Writes OUT/<read type>_<variant>.tsv.gz (training rows with p_species and p_samples, test rows with p) and appends
OUT/knobs.tsv.

    python3 ablate.py --build local/v15 --out ablate --read-types pe,se,pb,ont -t 4
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
from model_features import feature_columns  # noqa: E402

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
SET = "normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity"
COMPLEXITY = ["sample_log_taxa", "sample_low_identity", "sample_identity"]
REF = ["su_rate_ref", "lu_rate_ref", "lsu_rate_ref"]
VARIANTS = ["v15", "no_complexity", "no_depth"]
KEEP = ["meta_sample", "meta_read_pairs", "taxon", "truth", "sample_log_fragments"] + COMPLEXITY


def columns_of(variant, columns):
    cols = feature_columns(columns, SET)
    drop = {"v15": [], "no_complexity": COMPLEXITY, "no_depth": ["sample_log_fragments"], "no_ref": REF}[variant]
    return [c for c in cols if c not in drop]


def model(seed):
    return HistGradientBoostingClassifier(max_iter=250, learning_rate=0.1, max_leaf_nodes=63, min_samples_leaf=20,
                                          l2_regularization=1.0, class_weight="balanced", early_stopping=False,
                                          random_state=seed)


def f1_at(y, p, t, w):
    call = p >= t
    tp, fp, fn = (w * (call & (y == 1))).sum(), (w * (call & (y == 0))).sum(), (w * (~call & (y == 1))).sum()
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def choose_knob(y, p, w):
    grid = np.round(np.arange(0.30, 0.851, 0.01), 2)
    f = np.array([f1_at(y, p, t, w) for t in grid])
    return float(grid[f.argmax()]) if f.max() >= f1_at(y, p, 0.5, w) + 0.002 else 0.5


def out_of_fold(X, y, w, splits, seed):
    p = np.full(len(y), np.nan)
    for tr, te in splits:
        p[te] = model(seed).fit(X[tr], y[tr], sample_weight=w[tr]).predict_proba(X[te])[:, 1]
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("-t", "--threads", type=int, default=4)
    opts = ap.parse_args()
    warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
    os.makedirs(opts.out, exist_ok=True)
    knobs = os.path.join(opts.out, "knobs.tsv")
    for rt in opts.read_types.split(","):
        train = trainer.load_table(os.path.join(opts.build, "training", TABLES[rt]))
        test = trainer.load_table(os.path.join(opts.build, "test", TABLES[rt]))
        y = train["truth"].to_numpy()
        sc = trainer.scenario_of(train)
        w = np.where(sc == "", 1.0, 0.25)
        species = list(GroupKFold(opts.folds).split(train, y, train["taxon"].astype(str)))
        samples = list(GroupKFold(opts.folds).split(train, y, train["meta_sample"].astype(str)))
        for variant in opts.variants.split(","):
            began = time.time()
            cols = columns_of(variant, train.columns)
            X = train[cols].to_numpy(dtype=np.float64)
            with threadpool_limits(limits=opts.threads, user_api="openmp"):
                p_species = out_of_fold(X, y, w, species, opts.seed)
                p_samples = out_of_fold(X, y, w, samples, opts.seed)
                final = model(opts.seed).fit(X, y, sample_weight=w)
                p_test = final.predict_proba(test[cols].to_numpy(dtype=np.float64))[:, 1]
            knob = choose_knob(y, p_species, w)
            a = train[KEEP].copy()
            a.insert(0, "part", "training")
            a["meta_scenario"] = sc
            a["p_species"], a["p_samples"], a["p"] = p_species, p_samples, np.nan
            b = test[KEEP].copy()
            b.insert(0, "part", "test")
            b["meta_scenario"] = trainer.scenario_of(test)
            b["p_species"], b["p_samples"], b["p"] = np.nan, np.nan, p_test
            pd.concat([a, b]).to_csv(os.path.join(opts.out, f"{rt}_{variant}.tsv.gz"), sep="\t", index=False,
                                     float_format="%.6g")
            row = pd.DataFrame([{"read type": rt, "variant": variant, "features": len(cols), "knob": knob,
                                 "seconds": round(time.time() - began)}])
            row.to_csv(knobs, sep="\t", index=False, mode="a", header=not os.path.exists(knobs))
            print(f"{rt} {variant}: {len(cols)} features, knob {knob}, {time.time() - began:.0f} s", flush=True)


if __name__ == "__main__":
    main()
