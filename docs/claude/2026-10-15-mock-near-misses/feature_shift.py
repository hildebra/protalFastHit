#!/usr/bin/env python3
"""Where real samples' present taxa lie among the training's present taxa, feature by feature.

For each model feature: the training table's present rows (truth 1) of comparable evidence (fragments within
[--min-fragments, --max-fragments]) give the reference distribution; each real taxon of the dumps (protal
--profile_truth, <profile>.truth_annotated) with truth 1 and fragments in the same range gets its value's percentile
there. Prints per feature the real taxa's median percentile, the share below the 1st and above the 99th percentile,
and the training and real medians, sorted by how far the real taxa sit from the middle; then per real taxon its most
extreme features. Needs numpy and pandas.

usage: feature_shift.py MODEL.xml TRAINING_TABLE DUMP [DUMP ...] [--min-fragments 200] [--max-fragments 20000]
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd


def model_features(path):
    """The model's input fields: the PMML's MiningField names that are not the target."""
    text = open(path).read()
    names = re.findall(r'<MiningField name="([^"]+)"(?: usageType="([a-z]+)")?', text)
    return list(dict.fromkeys(n for n, use in names if use not in ("target", "predicted")))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model")
    ap.add_argument("training")
    ap.add_argument("dumps", nargs="+")
    ap.add_argument("--min-fragments", type=float, default=200)
    ap.add_argument("--max-fragments", type=float, default=20000)
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()

    features = model_features(args.model)
    train = pd.read_csv(args.training, sep="\t", low_memory=False)
    label = "truth" if "truth" in train.columns else "label"
    ref = train[(train[label] == 1) & train["fragments"].between(args.min_fragments, args.max_fragments)]
    features = [f for f in features if f in ref.columns]
    print(f"model features {len(features)}; training present rows in range {len(ref)} of {int((train[label] == 1).sum())}")

    real = []
    for dump in args.dumps:
        table = pd.read_csv(dump, sep="\t", low_memory=False)
        sample = os.path.basename(dump).split(".profile")[0]
        rows = table[(table["truth"] == 1) & table["fragments"].between(args.min_fragments, args.max_fragments)].copy()
        rows["sample"] = sample
        real.append(rows)
    real = pd.concat(real, ignore_index=True)
    print(f"real present taxa in range {len(real)}; called {int(real['prediction'].sum())}; "
          f"median probability {real['probability'].median():.3f}")

    pct = pd.DataFrame(index=real.index)
    out = []
    for f in features:
        ref_values = np.sort(ref[f].astype(float).to_numpy())
        values = real[f].astype(float).to_numpy()
        p = np.searchsorted(ref_values, values, side="left") / len(ref_values)
        p_right = np.searchsorted(ref_values, values, side="right") / len(ref_values)
        p = (p + p_right) / 2
        pct[f] = p
        out.append({"feature": f, "real_median_pct": float(np.median(p)), "below_p1": float((p < 0.01).mean()),
                    "above_p99": float((p > 0.99).mean()), "train_median": float(np.median(ref_values)),
                    "real_median": float(np.median(values)),
                    "train_p1": float(np.quantile(ref_values, 0.01)), "train_p99": float(np.quantile(ref_values, 0.99))})
    out = pd.DataFrame(out)
    out["off_centre"] = (out["real_median_pct"] - 0.5).abs() + out["below_p1"] + out["above_p99"]
    out = out.sort_values("off_centre", ascending=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(out.head(args.top).to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print("\nper taxon: the five features farthest out (percentile)")
    for i, row in real.iterrows():
        extreme = (pct.loc[i] - 0.5).abs().sort_values(ascending=False).head(5)
        parts = ", ".join(f"{f} {pct.loc[i, f]:.3f}" for f in extreme.index)
        print(f"{row['sample']}\t{row['taxon_name']}\tfragments {int(row['fragments'])}\tP {row['probability']:.3f}\t{parts}")


if __name__ == "__main__":
    sys.exit(main())
