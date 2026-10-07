#!/usr/bin/env python3
"""Could the knob follow a sample's complexity instead of being one for every sample? (r226 v14 predictions, no refit)

The global knob (pe 0.82, se 0.73) is chosen on training rows of which the soil scenarios are a large share; the
design's test set would rather have ~0.58. A sample's taxa with reads (rows of the profile; known to protal before
calling) tell a simple community from a soil one. Per bin of log10(taxa with reads), the knob of the highest F1 with
species held out on the training rows in that bin (rows weighted as in the fit: scenario rows 0.25; kept only if it
gains 0.002 over the global knob), applied to the design's test set and the scenarios' hold-out samples (by the final
model), against the global knob and 0.5.

    python3 knob_by_taxa.py --build local/v14 --read-types pe,se
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

SUFFIX = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
TABLE = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
         "ont": "training_data_ont.tsv"}
EDGES = [0, 2.5, 3.0, 3.5, 4.0, 9]


def table(path):
    df = pd.read_csv(path, sep="\t", usecols=["meta_sample", "meta_scenario", "taxon", "truth"], dtype={"meta_scenario": str})
    df["meta_scenario"] = df["meta_scenario"].fillna("")
    df["log_taxa"] = np.log10(df["meta_sample"].map(df["meta_sample"].value_counts()).astype(float))
    df["bin"] = np.digitize(df["log_taxa"], EDGES[1:-1])
    return df


def f1(y, call, w=None):
    w = np.ones(len(y)) if w is None else w
    tp, fp, fn = (w * (call & (y == 1))).sum(), (w * (call & (y == 0))).sum(), (w * (~call & (y == 1))).sum()
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", required=True)
    ap.add_argument("--read-types", default="pe,se")
    opts = ap.parse_args()
    grid = np.round(np.arange(0.30, 0.901, 0.01), 2)
    for rt in opts.read_types.split(","):
        ml = os.path.join(opts.build, "model_logs")
        with open(os.path.join(ml, f"trained_model{SUFFIX[rt]}.metrics.json")) as fh:
            knob = float(((json.load(fh).get("depth_knobs") or {}).get("global_knob")) or 0.5)
        train = table(os.path.join(opts.build, "training", TABLE[rt]))
        train["p"] = pd.read_csv(os.path.join(ml, f"trained_model{SUFFIX[rt]}.predictions.tsv.gz"), sep="\t",
                                 usecols=["p_species"])["p_species"].to_numpy()
        test = table(os.path.join(opts.build, "test", TABLE[rt]))
        p = np.full(len(test), np.nan)
        d = (test["meta_scenario"] == "").to_numpy()
        p[d] = pd.read_csv(os.path.join(ml, f"trained_model{SUFFIX[rt]}.test_predictions.tsv.gz"), sep="\t", usecols=["p"])["p"]
        p[~d] = pd.read_csv(os.path.join(ml, f"trained_model{SUFFIX[rt]}.scenario_predictions.tsv.gz"), sep="\t",
                            usecols=["p"])["p"]
        test["p"] = p
        w = np.where(train["meta_scenario"] == "", 1.0, 0.25)
        y = train["truth"].to_numpy()
        print(f"== {rt}: global knob {knob}")
        print("bin (log10 taxa)  training samples (design/scenario)  rows  knob  F1 there at the global knob -> at its own")
        knobs = {}
        for b in range(len(EDGES) - 1):
            m = (train["bin"] == b).to_numpy()
            if not m.any():
                continue
            f = np.array([f1(y[m], train["p"].to_numpy()[m] >= t, w[m]) for t in grid])
            base = f1(y[m], train["p"].to_numpy()[m] >= knob, w[m])
            k = float(grid[f.argmax()]) if f.max() >= base + 0.002 else knob
            knobs[b] = k
            sm = train[m].groupby("meta_sample")["meta_scenario"].first()
            print(f"  {EDGES[b]:.1f}-{EDGES[b + 1]:.1f}  {int((sm == '').sum())}/{int((sm != '').sum())}  {int(m.sum())}  {k:.2f}  "
                  f"{base:.4f} -> {f1(y[m], train['p'].to_numpy()[m] >= k, w[m]):.4f}")
        k_row = test["bin"].map(lambda b: knobs.get(b, knob)).to_numpy()
        print("set  samples  F1 at 0.5 / at the global knob / at the knob of its bin (FP, FN)")
        for label, g in [("design test", test[d])] + [(f"{s} hold-out", x) for s, x in test[~d].groupby("meta_scenario")]:
            yy, pp, kk = g["truth"].to_numpy(), g["p"].to_numpy(), k_row[g.index]
            call = pp >= kk
            print(f"  {label}  {g['meta_sample'].nunique()}  {f1(yy, pp >= 0.5):.4f} / {f1(yy, pp >= knob):.4f} / "
                  f"{f1(yy, call):.4f} ({int((call & (yy == 0)).sum())}, {int((~call & (yy == 1)).sum())})")


if __name__ == "__main__":
    main()
