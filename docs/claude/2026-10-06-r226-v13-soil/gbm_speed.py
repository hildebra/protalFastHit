#!/usr/bin/env python3
"""How fast gradient boosting fits the r226 v13 tables, and what cheaper settings cost: one fit on four fifths of a
read type's training rows (one species fold) at several OpenMP thread counts, and the F1 of settings with fewer rounds
(at a higher learning rate) or fewer histogram bins, with species held out (5 folds) and on the build's test table
(the design's test set and the soil hold-out samples), the scenarios' rows at weight 0.25 and the sample's depth
rounded to 0.25 as varied scenario depths will make it.

    python3 gbm_speed.py --training local/v13/training --test local/v13/test --read-type pe --out speed_pe.tsv
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import machine_learning_cmdline as trainer  # noqa: E402
from model_features import DEFAULT_FEATURE_SET, feature_columns  # noqa: E402

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
BASE = dict(max_iter=500, learning_rate=0.05, max_leaf_nodes=63, min_samples_leaf=20, l2_regularization=1.0,
            class_weight="balanced", early_stopping=False, random_state=1)
SETTINGS = {"500 rounds, rate 0.05 (default)": {}, "250 rounds, rate 0.1": {"max_iter": 250, "learning_rate": 0.1},
            "150 rounds, rate 0.15": {"max_iter": 150, "learning_rate": 0.15},
            "500 rounds, 63 bins": {"max_bins": 63}}


def f1(y, p, knob=0.5):
    call = p >= knob
    tp, fp, fn = (call & (y == 1)).sum(), (call & (y == 0)).sum(), (~call & (y == 1)).sum()
    return 2 * tp / (2 * tp + fp + fn)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--training", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--read-type", default="pe")
    ap.add_argument("--threads", default="1,2,3,6")
    ap.add_argument("--out", required=True)
    opts = ap.parse_args()
    train = trainer.load_table(os.path.join(opts.training, TABLES[opts.read_type]))
    test = trainer.load_table(os.path.join(opts.test, TABLES[opts.read_type]))
    cols = feature_columns(train.columns, DEFAULT_FEATURE_SET)
    X, Xt = train[cols].to_numpy(dtype=np.float64), test[cols].to_numpy(dtype=np.float64)
    j = cols.index("sample_log_fragments")
    X[:, j], Xt[:, j] = np.round(X[:, j] * 4) / 4, np.round(Xt[:, j] * 4) / 4
    y, yt = train["truth"].to_numpy(), test["truth"].to_numpy()
    w = np.where(trainer.scenario_of(train) != "", 0.25, 1.0)
    sct = trainer.scenario_of(test)
    splits = list(GroupKFold(5).split(X, y, train["taxon"].astype(str)))
    rows = []
    tr = splits[0][0]
    for threads in (int(t) for t in opts.threads.split(",")):
        with threadpool_limits(limits=threads, user_api="openmp"):
            t0 = time.time()
            HistGradientBoostingClassifier(**BASE).fit(X[tr], y[tr], sample_weight=w[tr])
            rows.append({"what": f"one fit, {threads} threads", "seconds": round(time.time() - t0, 1)})
            print(rows[-1], flush=True)
    threads = max(int(t) for t in opts.threads.split(","))
    for name, change in SETTINGS.items():
        params = {**BASE, **change}
        p = np.full(len(y), np.nan)
        t0 = time.time()
        with threadpool_limits(limits=threads, user_api="openmp"):
            for a, b in splits:
                p[b] = HistGradientBoostingClassifier(**params).fit(X[a], y[a], sample_weight=w[a]).predict_proba(X[b])[:, 1]
            seconds = (time.time() - t0) / len(splits)
            pt = HistGradientBoostingClassifier(**params).fit(X, y, sample_weight=w).predict_proba(Xt)[:, 1]
        row = {"what": name, "seconds": round(seconds, 1), "F1 species held out": f1(y, p),
               "F1 design test": f1(yt[sct == ""], pt[sct == ""])}
        for sc in ("soil", "soil_shallow"):
            if (sct == sc).any():
                row[f"F1 {sc} hold-out"] = f1(yt[sct == sc], pt[sct == sc])
        rows.append(row)
        print(row, flush=True)
    pd.DataFrame(rows).to_csv(opts.out, sep="\t", index=False, float_format="%.4f")


if __name__ == "__main__":
    main()
