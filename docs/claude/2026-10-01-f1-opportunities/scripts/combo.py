#!/usr/bin/env python3
"""combo.py [READ_TYPES...] - do the two candidates with consistent small gains stack? The forest refitted with the base
features, + the quality-calibrated divergence (qual_calib.py), + the conservation pattern of the hit genes
(features_exp.py), and + both; scored as own_cluster.py. Writes results/combo.md."""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
argv = sys.argv[1:]
sys.argv = [sys.argv[0]]
import qual_calib as Q  # noqa: E402
import features_exp as F  # noqa: E402

O, X = Q.O, Q.X
READ_TYPES = argv or ["pe", "ont"]


def main():
    lines = ["# Stacking the quality-calibrated divergence and the conservation pattern (combo.py)", ""]
    for rt in READ_TYPES:
        train, test = Q.with_q("training", rt), Q.with_q("test", rt)
        cons = F.add_conserv(train, "training", rt)
        F.add_conserv(test, "test", rt)
        lines += [f"## {rt}", "", "| features | CV F1 | test F1 | test FP | test FN | delta (interval) |", "|---|---|---|---|---|---|"]
        base_calls = {}
        for label, feats in (("base", O.BASE), ("+ excess", O.BASE + Q.FEATURES), ("+ conservation pattern", O.BASE + cons),
                             ("+ both", O.BASE + Q.FEATURES + cons)):
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
    with open(os.path.join(HERE, "..", "results", "combo.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
