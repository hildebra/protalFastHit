#!/usr/bin/env python3
"""Which change made v14's soil better: its training samples, or its model? Refits on both builds' tables.

Each model (--models) is fitted on each build's training table (--builds NAME=DIR,...) and scored on the same rows
for all: the design's test set (the same 22,305 rows in v13 and v14), and every build's scenario hold-out samples
(v13: 2 per scenario at the preset's depth; v14: 3 at depths varied 0.5-2x). The knob is chosen as the trainer chooses
its global one: the threshold of the highest F1 with species held out (5 folds grouped by taxon) over the training
rows, weighted as in the fit (scenario rows 0.25), kept if it gains 0.002 over 0.5.

Models:
  forest    the v13 build's: normalized+adjacency+distance+depth+divergence+unfiltered, 64 trees, 512 leaves (pe, se)
            or 128 (pb, ont), max_features sqrt, balanced classes
  gbm+ref   the v14 build's: the same features + ref (su/lu/lsu_rate_ref), 500 rounds at 0.05, 63 leaves, min leaf 20,
            L2 1, balanced classes, no early stopping
  gbm       gbm+ref without the ref features
  gbm250    gbm+ref at the trainer's defaults since 15006c2 (250 rounds at 0.1)

Writes OUT/sets.tsv (each set pooled) and OUT/samples.tsv (each hold-out sample), appending as it goes.

    python3 cross_fit.py --builds v13=local/v13,v14=local/v14 --out cross -t 6
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.model_selection import GroupKFold
from threadpoolctl import threadpool_limits

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import machine_learning_cmdline as trainer  # noqa: E402
from model_features import feature_columns  # noqa: E402

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
BASE = "normalized+adjacency+distance+depth+divergence+unfiltered"
LEAVES = {"pe": 512, "se": 512, "pb": 128, "ont": 128}


def features(name, columns):
    return feature_columns(columns, BASE if name in ("forest", "gbm") else BASE + "+ref")


def model(name, rt, seed, threads):
    if name == "forest":
        return RandomForestClassifier(n_estimators=64, max_leaf_nodes=LEAVES[rt], min_samples_leaf=1, max_features="sqrt",
                                      class_weight="balanced", random_state=seed, n_jobs=threads)
    rounds, rate = (250, 0.1) if name == "gbm250" else (500, 0.05)
    return HistGradientBoostingClassifier(max_iter=rounds, learning_rate=rate, max_leaf_nodes=63, min_samples_leaf=20,
                                          l2_regularization=1.0, class_weight="balanced", early_stopping=False,
                                          random_state=seed)


def f1_at(y, p, t, w=None):
    w = np.ones(len(y)) if w is None else w
    call = p >= t
    tp, fp, fn = (w * (call & (y == 1))).sum(), (w * (call & (y == 0))).sum(), (w * (~call & (y == 1))).sum()
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def choose_knob(y, p, w):
    grid = np.round(np.arange(0.30, 0.851, 0.01), 2)
    f = np.array([f1_at(y, p, t, w) for t in grid])
    return float(grid[f.argmax()]) if f.max() >= f1_at(y, p, 0.5, w) + 0.002 else 0.5


def best(y, p):
    grid = np.round(np.arange(0.05, 0.96, 0.01), 2)
    f = np.array([f1_at(y, p, t) for t in grid])
    return float(f.max()), float(grid[f.argmax()])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--builds", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    ap.add_argument("--models", default="forest,gbm+ref,gbm")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("-t", "--threads", type=int, default=6)
    opts = ap.parse_args()
    warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
    os.makedirs(opts.out, exist_ok=True)
    builds = [s.split("=", 1) for s in opts.builds.split(",")]
    sets_out, samples_out = [], []
    for rt in opts.read_types.split(","):
        train = {b: trainer.load_table(os.path.join(d, "training", TABLES[rt])) for b, d in builds}
        test = {b: trainer.load_table(os.path.join(d, "test", TABLES[rt])) for b, d in builds}
        # the scoring sets: the design's test rows once (identical in the builds), each build's hold-out samples
        first = builds[0][0]
        evals = [("design test", test[first][trainer.scenario_of(test[first]) == ""])]
        for b, _ in builds:
            sc = trainer.scenario_of(test[b])
            for scenario in sorted(set(sc) - {""}):
                evals.append((f"{scenario} hold-out {b}", test[b][sc == scenario]))
        for tb, _ in builds:
            df = train[tb]
            y = df["truth"].to_numpy()
            sc = trainer.scenario_of(df)
            w = np.where(sc == "", 1.0, 0.25)
            splits = list(GroupKFold(opts.folds).split(df, y, df["taxon"].astype(str)))
            for name in opts.models.split(","):
                began = time.time()
                cols = features(name, df.columns)
                X = df[cols].to_numpy(dtype=np.float64)
                p = np.full(len(y), np.nan)
                with threadpool_limits(limits=opts.threads, user_api="openmp"):
                    for tr, te in splits:
                        p[te] = model(name, rt, opts.seed, opts.threads).fit(X[tr], y[tr], sample_weight=w[tr]) \
                            .predict_proba(X[te])[:, 1]
                    final = model(name, rt, opts.seed, opts.threads).fit(X, y, sample_weight=w)
                    knob = choose_knob(y, p, w)
                    scored = [(label, rows["truth"].to_numpy(), final.predict_proba(rows[cols].to_numpy(dtype=np.float64))[:, 1],
                               rows) for label, rows in evals]
                for scenario in sorted(set(sc) - {""}):
                    m = sc == scenario
                    scored.append((f"{scenario} hold-in {tb}, species held out", y[m], p[m], df[m]))
                scored.append((f"design hold-in {tb}, species held out", y[sc == ""], p[sc == ""], df[sc == ""]))
                for label, yy, pp, rows in scored:
                    call = pp >= knob
                    bf, bt = best(yy, pp)
                    sets_out.append({"read type": rt, "trained on": tb, "model": name, "knob": knob, "set": label,
                                     "taxa": len(yy), "FP": int((call & (yy == 0)).sum()), "FN": int((~call & (yy == 1)).sum()),
                                     "F1": f1_at(yy, pp, knob), "F1 at 0.5": f1_at(yy, pp, 0.5), "best F1": bf, "at": bt})
                    if "hold-out" not in label:
                        continue
                    for sample, idx in rows.groupby("meta_sample").indices.items():
                        ys, ps = yy[idx], pp[idx]
                        cs = ps >= knob
                        samples_out.append({"read type": rt, "trained on": tb, "model": name, "knob": knob, "set": label,
                                            "sample": sample,
                                            "log10 fragments": float(rows["sample_log_fragments"].iloc[idx[0]]),
                                            "FP": int((cs & (ys == 0)).sum()), "FN": int((~cs & (ys == 1)).sum()),
                                            "F1": f1_at(ys, ps, knob), "F1 at 0.5": f1_at(ys, ps, 0.5)})
                print(f"{rt} trained on {tb}, {name}: knob {knob}, {time.time() - began:.0f} s; " +
                      ", ".join(f"{r['set']} {r['F1']:.4f}" for r in sets_out
                                if r["read type"] == rt and r["trained on"] == tb and r["model"] == name), flush=True)
                pd.DataFrame(sets_out).to_csv(os.path.join(opts.out, "sets.tsv"), sep="\t", index=False, float_format="%.5g")
                pd.DataFrame(samples_out).to_csv(os.path.join(opts.out, "samples.tsv"), sep="\t", index=False,
                                                 float_format="%.5g")


if __name__ == "__main__":
    main()
