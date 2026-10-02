#!/usr/bin/env python3
"""Where the r226 build's models err: false positives and negatives of each read type (species held out at
knob 0.5, and the independent test set), by their source, depth, fragments and identity; and a threshold
per sample depth, chosen on the training rows (species held out) and applied to the test set, whose depths
differ (interpolated in log10 depth), as an out-of-sample check. Usage: errors.py MODEL_LOGS_DIR"""
import sys

import numpy as np
import pandas as pd

logs = sys.argv[1]
pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 30)
BINS, LABELS = [0, 1, 2, 5, 10, 100, 1e12], ["1", "2", "3-5", "6-10", "11-100", ">100"]
THRESHOLDS = np.round(np.arange(0.10, 0.951, 0.02), 2)


def f1(truth, call):
    tp, fp, fn = int((truth & call).sum()), int((~truth & call).sum()), int((truth & ~call).sum())
    return 2 * tp / max(1, 2 * tp + fp + fn), tp, fp, fn


def load(suffix, name):
    df = pd.read_csv(f"{logs}/trained_model{suffix}.{name}.tsv.gz", sep="\t", low_memory=False)
    df["is_present"] = df["truth"].astype(str).str.lower().isin(("1", "true"))
    novel = df["meta_novel_level"].fillna("").astype(str).replace("nan", "")
    df["source"] = np.where(novel != "", "near held-out " + novel, "near a species the db has")
    return df


def report(df, pcol, label):
    truth, call = df["is_present"], df[pcol] >= 0.5
    print(f"\n=================== {label} ({len(df)} taxa, {df['meta_sample'].nunique()} samples)")
    print("F1 %.4f  TP %d FP %d FN %d" % f1(truth, call))
    fp, fn = df[~truth & call], df[truth & ~call]
    print("\nFP by source x deepest rank shared with a simulated species:")
    print(pd.crosstab(fp["source"], fp["meta_relative_rank"].fillna("-"), margins=True))
    if "fragments" in df:
        print("\nFP by source x fragments:")
        print(pd.crosstab(fp["source"], pd.cut(fp["fragments"], BINS, labels=LABELS), margins=True))
        perfect = fp["top_identity"] >= 0.9999
        near_db = fp["source"].eq("near a species the db has")
        print(f"FP with a read identical to the taxon's gene (top_identity 1): {perfect.sum()} of {len(fp)}; "
              f"near a species the db has: {(perfect & near_db).sum()} of {near_db.sum()}")
        print("\nFN by simulated from the representative (1) or another genome (0) x fragments:")
        print(pd.crosstab(fn["meta_rep_genome"].fillna("-").astype(str), pd.cut(fn["fragments"], BINS, labels=LABELS),
                          margins=True))
        strains = df[truth & df["meta_rep_genome"].eq(0)]
        print("present taxa simulated from another genome, identity quantiles 5/25/50%%: found %s | missed %s" % (
            np.round(strains.loc[strains[pcol] >= 0.5, "identity"].quantile([.05, .25, .5]).to_numpy(), 4),
            np.round(strains.loc[strains[pcol] < 0.5, "identity"].quantile([.05, .25, .5]).to_numpy(), 4)))
        many = fn[fn["fragments"] > 10]
        if len(many):
            print(f"FN with more than 10 fragments: {len(many)}, from another genome {int((many['meta_rep_genome'] == 0).sum())}, "
                  f"median identity {many['identity'].median():.4f}, median top_identity {many['top_identity'].median():.4f}")


def per_depth(df, pcol):
    """Rows (depth, samples, FP per sample at 0.5, F1 at 0.5, best threshold, F1 at it)."""
    rows = []
    for depth, g in df.groupby("meta_read_pairs"):
        truth = g["is_present"]
        base = f1(truth, g[pcol] >= 0.5)
        best = max((f1(truth, g[pcol] >= th)[0], th) for th in THRESHOLDS)
        rows.append((int(depth), g["meta_sample"].nunique(), base[2] / g["meta_sample"].nunique(), base[0], best[1], best[0]))
    return pd.DataFrame(rows, columns=["depth", "samples", "FP_per_sample", "F1_0.5", "best_threshold", "F1_best"])


for t, suffix in (("pe", ""), ("se", "_se"), ("pb", "_pb"), ("ont", "_ont")):
    train, test = load(suffix, "predictions"), load(suffix, "test_predictions")
    report(train, "p_species", f"{t}, species held out")
    report(test, "p", f"{t}, independent test set")
    knobs = per_depth(train, "p_species")
    print(f"\n{t}, species held out, by depth:")
    print(knobs.to_string(index=False))
    print(f"\n{t}, test set, by depth:")
    print(per_depth(test, "p").to_string(index=False))
    # Out of sample: the training rows' best threshold per depth, interpolated in log10 depth, on the test set.
    x, y = np.log10(knobs["depth"].to_numpy(float)), knobs["best_threshold"].to_numpy(float)
    th = np.interp(np.log10(test["meta_read_pairs"].to_numpy(float)), x, y)
    base, tuned = f1(test["is_present"], test["p"] >= 0.5), f1(test["is_present"], test["p"] >= th)
    print(f"\n{t}, test set with the training's threshold per depth (interpolated): F1 {tuned[0]:.4f} "
          f"(TP {tuned[1]}, FP {tuned[2]}, FN {tuned[3]}) against {base[0]:.4f} (TP {base[1]}, FP {base[2]}, FN {base[3]}) at 0.5")
