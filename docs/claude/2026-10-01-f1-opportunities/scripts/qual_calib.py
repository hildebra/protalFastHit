#!/usr/bin/env python3
"""qual_calib.py [READ_TYPES...] - divergence beyond what the base qualities explain, as a feature.

A read's differences from the reference are sequencing errors plus the divergence of its genome from the reference
(a strain's, or a relative's). The base qualities predict the errors: a read's excess = its differences per aligned
base (X + I + D over M + X + I + D) less the mean error probability of its bases (10^(-Q/10)). Per (sample, taxon),
from the SAMs of 0.7.1's pipeline: excess_median (the median over its primary records) and excess_high_share (the
share of records whose excess is above 0.02). Long reads, whose error rate varies from read to read (pbsim3's
quality model), should gain most. The forest refitted with the base features and with these added, as
own_cluster.py does. Writes results/qual_calib.md.
"""
import collections
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.argv = [sys.argv[0]] + sys.argv[1:]
import own_cluster as O  # noqa: E402  (paths, SAM lookup, the forest and scores)

X = O.X
READ_TYPES = sys.argv[1:] or ["pe", "ont"]
FEATURES = ["excess_median", "excess_high_share"]
CACHE = os.path.join(os.environ["FPEXP"], "qual")
ERR = np.array([10 ** (-(q - 33) / 10) for q in range(128)])


def features(split, rt):
    import pandas as pd
    path = os.path.join(CACHE, f"{split}_{rt}.tsv")
    if os.path.exists(path):
        return pd.read_csv(path, sep="\t")
    rows = []
    for sample in sorted(X.load(split, rt)["meta_sample"].unique()):
        sam = O.sam_of(split, rt, sample)
        if not os.path.exists(sam):
            print(f"no SAM for {sample}", file=sys.stderr)
            continue
        per = collections.defaultdict(list)
        for f in X.sam_records(sam):
            flag = int(f[1])
            if flag & (4 | 256 | 2048) or len(f) < 11 or f[10] == "*":
                continue
            c = X.cigar_counts(f[5])
            aligned = c["M"] + c["X"] + c["I"] + c["D"]
            if aligned == 0:
                continue
            q = np.frombuffer(f[10].encode(), dtype=np.uint8)
            expected = float(ERR[np.clip(q, 0, 127)].mean())
            per[f[2].split("_")[0]].append((c["X"] + c["I"] + c["D"]) / aligned - expected)
        for taxid, ex in per.items():
            ex = np.array(ex)
            rows.append({"meta_sample": sample, "taxon": int(taxid), "excess_median": float(np.median(ex)),
                         "excess_high_share": float((ex > 0.02).mean())})
    df = pd.DataFrame(rows)
    os.makedirs(CACHE, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    return df


def with_q(split, rt):
    t = X.load(split, rt).merge(features(split, rt), on=["meta_sample", "taxon"], how="left")
    t["excess_median"] = t["excess_median"].fillna(t["excess_median"].median())
    t["excess_high_share"] = t["excess_high_share"].fillna(0.0)
    return t


def main():
    lines = ["# Divergence beyond the base qualities (qual_calib.py)", "",
             f"Seeds {O.SEEDS}; as own_cluster.md.", ""]
    for rt in READ_TYPES:
        train, test = with_q("training", rt), with_q("test", rt)
        pres, absent = test[test["truth"] == 1], test[test["truth"] == 0]
        lines += [f"## {rt}", "", f"excess_median: present taxa {pres['excess_median'].median():.4f} (with a missing "
                  f"congener in the sample {pres.loc[pres['meta_novel_congener'] == 1, 'excess_median'].median():.4f}), "
                  f"absent taxa {absent['excess_median'].median():.4f}", "",
                  "| features | CV F1 | test F1 | test FP | test FN | delta (interval) |", "|---|---|---|---|---|---|"]
        base_calls = {}
        for label, feats in (("base", O.BASE), ("base + excess", O.BASE + FEATURES)):
            cv, tf, fp, fn, diffs = [], [], [], [], []
            for seed in O.SEEDS:
                p_cv = X.cv_predict(train, feats, seed)
                call = X.test_predict(train, test, feats, seed) >= X.KNOB
                y = test["truth"].to_numpy()
                cv.append(X.scores(train["truth"].to_numpy(), p_cv >= X.KNOB)["F1"])
                s = X.scores(y, call)
                tf.append(s["F1"]); fp.append(s["FP"]); fn.append(s["FN"])
                if label == "base":
                    base_calls[seed] = call
                else:
                    diffs.append(X.bootstrap_diff(test, base_calls[seed], call, n=500, seed=seed))
            d = np.array(diffs).mean(axis=0) if diffs else (0, 0, 0)
            lines.append(f"| {label} | {np.mean(cv):.4f} | {np.mean(tf):.4f} | {np.mean(fp):.1f} | {np.mean(fn):.1f} | "
                         f"{d[0]:+.4f} ({d[1]:+.4f}, {d[2]:+.4f}) |")
            print(lines[-1], flush=True)
        lines.append("")
    with open(os.path.join(HERE, "..", "results", "qual_calib.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
