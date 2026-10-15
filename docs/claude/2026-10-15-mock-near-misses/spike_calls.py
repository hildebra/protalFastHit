#!/usr/bin/env python3
"""The calls of one spike mode against another (two of spike_compare.sh's folders of protal training
dumps, <sample>.profile.truth_annotated): per sample TP/FP/FN under each, the taxa whose call changes (with their
probabilities), and the totals over the mock samples and over the simulated samples (names starting with rl, pb_ or
ont_ are simulated). Needs pandas.

usage: spike_calls.py FOLDER BASELINE_FOLDER   (e.g. auto off: the first is "rule", the baseline "kept")
"""
import glob
import os
import sys

import pandas as pd


def calls(path):
    t = pd.read_csv(path, sep="\t", usecols=["truth", "prediction", "probability", "taxon", "taxon_name"], low_memory=False)
    return t.set_index("taxon")


def counts(t):
    tp = int(((t.truth == 1) & (t.prediction == 1)).sum())
    fp = int(((t.truth == 0) & (t.prediction == 1)).sum())
    fn = int(((t.truth == 1) & (t.prediction == 0)).sum())
    return tp, fp, fn


def f1(tp, fp, fn):
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def main():
    rule, keep = sys.argv[1], sys.argv[2]
    totals = {"mock": [0] * 6, "sim": [0] * 6}
    changed_samples = 0
    print("sample\tTP/FP/FN kept\tTP/FP/FN rule\tchanges")
    for path in sorted(glob.glob(os.path.join(keep, "*.profile.truth_annotated"))):
        sample = os.path.basename(path).split(".profile")[0]
        other = os.path.join(rule, os.path.basename(path))
        if not os.path.exists(other):
            print(f"{sample}\tmissing in {rule}")
            continue
        a, b = calls(path), calls(other)
        ca, cb = counts(a), counts(b)
        kind = "sim" if sample.startswith(("rl", "pb_", "ont_")) else "mock"
        for i, v in enumerate(ca + cb):
            totals[kind][i] += v
        joined = a.join(b, how="outer", lsuffix="_keep", rsuffix="_rule")
        flips = joined[(joined.prediction_keep.fillna(0) != joined.prediction_rule.fillna(0))]
        moved = joined[(joined.probability_keep.fillna(0) - joined.probability_rule.fillna(0)).abs() > 1e-9]
        parts = [f"{r.taxon_name_keep if isinstance(r.taxon_name_keep, str) else r.taxon_name_rule} "
                 f"(truth {int(r.truth_keep if r.truth_keep == r.truth_keep else r.truth_rule)}) "
                 f"P {r.probability_keep:.3f} -> {r.probability_rule:.3f}" for _, r in flips.iterrows()]
        if len(moved):
            changed_samples += 1
        if parts or kind == "mock":
            print(f"{sample}\t{ca[0]}/{ca[1]}/{ca[2]}\t{cb[0]}/{cb[1]}/{cb[2]}\t"
                  f"{len(moved)} taxa moved; calls changed: {'; '.join(parts) if parts else 'none'}")
    for kind, v in totals.items():
        print(f"# {kind}: kept TP/FP/FN {v[0]}/{v[1]}/{v[2]} F1 {f1(*v[:3]):.4f}; rule {v[3]}/{v[4]}/{v[5]} F1 {f1(*v[3:]):.4f}")
    print(f"# samples with any probability moved: {changed_samples}")


if __name__ == "__main__":
    sys.exit(main())
