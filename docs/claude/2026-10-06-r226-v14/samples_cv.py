#!/usr/bin/env python3
"""Does a scenario sample's depth still stand for its scenario? Whole samples held out, with and without the depth and
with the sample's complexity as features (r226 v14 tables, gradient boosting as v14 trained it).

The trainer's "samples held out" estimate (GroupKFold over meta_sample, 5 folds; p_samples) puts pe's shallowest
shallow-soil sample at F1 0.758 against 0.920 with species held out: held out, its depth (4.86) is below every other
soil sample's, and the model scores it like the design's samples of that depth. Each variant here is scored the same
way (the same folds), per scenario sample, and its final model on both builds' hold-out samples:

  v14         the build's model: normalized+adjacency+distance+depth+divergence+unfiltered+ref, 500 rounds at 0.05,
              63 leaves (reproduces the build's p_samples)
  no depth    without sample_log_fragments
  +sample     + the sample's complexity, known before calling: sample_log_taxa (log10 taxa with reads),
              sample_low_identity (the sample's fragment share on low-identity bases), sample_identity (its taxa's
              fragment-weighted median identity) (soil_experiments.add_context, r226 v13 report)
  +sample no depth   both

Calls at the build's knob (--knob, 0.82 for pe), at 0.5, and at the best threshold of each set's rows (the variants'
own knobs would need the species folds too: 5 more fits each).

    python3 samples_cv.py --build local/v14 --other local/v13 --read-types pe -t 6
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
sys.path.insert(0, os.path.join(HERE, "..", "2026-10-06-r226-v13-soil"))
import machine_learning_cmdline as trainer  # noqa: E402
from model_features import feature_columns  # noqa: E402
from soil_experiments import add_context  # noqa: E402

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
SET = "normalized+adjacency+distance+depth+divergence+unfiltered+ref"
SAMPLE = ["sample_log_taxa", "sample_low_identity", "sample_identity"]
VARIANTS = ["v14", "no depth", "+sample", "+sample no depth"]


def columns(variant, cols):
    out = feature_columns(cols, SET)
    if "no depth" in variant:
        out = [c for c in out if c != "sample_log_fragments"]
    if "+sample" in variant:
        out += SAMPLE
    return out


def model(seed):
    return HistGradientBoostingClassifier(max_iter=500, learning_rate=0.05, max_leaf_nodes=63, min_samples_leaf=20,
                                          l2_regularization=1.0, class_weight="balanced", early_stopping=False,
                                          random_state=seed)


def f1(y, c):
    tp, fp, fn = int((c & (y == 1)).sum()), int((c & (y == 0)).sum()), int((~c & (y == 1)).sum())
    return (2 * tp / (2 * tp + fp + fn) if tp else 0.0), fp, fn


def best(y, p):
    grid = np.round(np.arange(0.05, 0.96, 0.01), 2)
    f = [f1(y, p >= t)[0] for t in grid]
    i = int(np.argmax(f))
    return f[i], float(grid[i])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", required=True)
    ap.add_argument("--other", help="another build whose hold-out samples the final models score too")
    ap.add_argument("--read-types", default="pe")
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--knob", type=float, default=0.82)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("-t", "--threads", type=int, default=6)
    opts = ap.parse_args()
    warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
    os.makedirs(opts.out, exist_ok=True)
    out = []
    for rt in opts.read_types.split(","):
        df = add_context(trainer.load_table(os.path.join(opts.build, "training", TABLES[rt])))
        tests = {"v14 hold-out": add_context(trainer.load_table(os.path.join(opts.build, "test", TABLES[rt])))}
        if opts.other:
            tests["v13 hold-out"] = add_context(trainer.load_table(os.path.join(opts.other, "test", TABLES[rt])))
        y = df["truth"].to_numpy()
        sc = trainer.scenario_of(df)
        w = np.where(sc == "", 1.0, 0.25)
        splits = list(GroupKFold(5).split(df, y, df["meta_sample"].astype(str).to_numpy()))
        for variant in opts.variants.split(","):
            began = time.time()
            cols = columns(variant, list(df.columns))
            X = df[cols].to_numpy(dtype=np.float64)
            p = np.full(len(y), np.nan)
            with threadpool_limits(limits=opts.threads, user_api="openmp"):
                for tr, te in splits:
                    p[te] = model(opts.seed).fit(X[tr], y[tr], sample_weight=w[tr]).predict_proba(X[te])[:, 1]
                final = model(opts.seed).fit(X, y, sample_weight=w)
                scored = {k: (t, final.predict_proba(t[cols].to_numpy(dtype=np.float64))[:, 1]) for k, t in tests.items()}
            sets = [("design, samples held out", df[sc == ""], p[sc == ""])]
            for s in sorted(set(sc) - {""}):
                sets.append((f"{s} hold-in, samples held out", df[sc == s], p[sc == s]))
            for k, (t, pt) in scored.items():
                tsc = trainer.scenario_of(t)
                if k == "v14 hold-out":
                    sets.append(("design test", t[tsc == ""], pt[tsc == ""]))
                for s in sorted(set(tsc) - {""}):
                    sets.append((f"{s} {k}", t[tsc == s], pt[tsc == s]))
            for label, rows, pp in sets:
                yy = rows["truth"].to_numpy()
                a, fp, fn = f1(yy, pp >= opts.knob)
                bf, bt = best(yy, pp)
                out.append({"read type": rt, "variant": variant, "set": label, "sample": "all", "log10 fragments": np.nan,
                            "F1": a, "FP": fp, "FN": fn, "F1 at 0.5": f1(yy, pp >= 0.5)[0], "best F1": bf, "at": bt})
                if label.startswith("design"):
                    continue
                for sample, idx in rows.groupby("meta_sample").indices.items():
                    ys, ps = yy[idx], pp[idx]
                    a, fp, fn = f1(ys, ps >= opts.knob)
                    out.append({"read type": rt, "variant": variant, "set": label, "sample": sample,
                                "log10 fragments": float(rows["sample_log_fragments"].iloc[idx[0]]), "F1": a, "FP": fp,
                                "FN": fn, "F1 at 0.5": f1(ys, ps >= 0.5)[0], "best F1": np.nan, "at": np.nan})
            print(f"{rt} {variant}: {time.time() - began:.0f} s; " + ", ".join(
                f"{r['set']} {r['F1']:.4f}" for r in out if r["variant"] == variant and r["read type"] == rt and r["sample"] == "all"),
                flush=True)
            pd.DataFrame(out).to_csv(os.path.join(opts.out, "samples_cv.tsv"), sep="\t", index=False, float_format="%.5g")


if __name__ == "__main__":
    main()
