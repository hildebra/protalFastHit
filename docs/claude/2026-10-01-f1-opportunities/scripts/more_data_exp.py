#!/usr/bin/env python3
"""more_data_exp.py [READ_TYPES...] - does twice the training data raise F1? 0.7.1's forest (the trainer's settings, the base
features) fitted on the pipeline's training table (A, seed 1), on a second collection of the same design (B,
more_training.sh: seed 101, ~/bench071/V071_training2) and on both, each scored on the pipeline's independent test
set; seeds 1-3, the test difference against A with a paired bootstrap over test samples. Writes
results/more_data_exp.md.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
V071 = os.path.expanduser(os.environ.get("V071", "~/bench071/V071"))
SECOND = os.path.expanduser(os.environ.get("SECOND", "~/bench071/V071_training2"))
os.environ["V2"] = V071
sys.path.insert(0, os.path.join(HERE, "..", "..", "2026-10-01-alignment-features"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts"))
import exp_lib as X  # noqa: E402
from model_features import NORMALIZED_FEATURES as BASE  # noqa: E402

READ_TYPES = sys.argv[1:] or ["pe", "se", "pb", "ont"]
SEEDS = [1, 2, 3]


def load_second(rt):
    df = pd.read_csv(os.path.join(SECOND, X.TABLE[rt]), sep="\t", float_precision="round_trip", low_memory=False)
    df["truth"] = df["truth"].astype(int)
    df["genus"] = df["taxon_name"].str.split(" ").str[0].str[3:]
    df["meta_sample"] = "B_" + df["meta_sample"].astype(str)  # distinct from A's sample names
    return df


def main():
    lines = ["# Twice the training data (more_data_exp.py)", "",
             "Test F1 on the pipeline's independent test set, mean of seeds 1-3; delta against A with a paired bootstrap "
             "over test samples (mean, 95% interval).", ""]
    for rt in READ_TYPES:
        a, b, test = X.load("training", rt), load_second(rt), X.load("test", rt)
        y = test["truth"].to_numpy()
        lines += [f"## {rt} (A {len(a)} taxa, B {len(b)}, test {len(test)})", "",
                  "| training | test F1 | FP | FN | delta against A (interval) |", "|---|---|---|---|---|"]
        calls_a = {}
        for label, train in (("A (the pipeline's)", a), ("B (seed 101)", b), ("A + B", pd.concat([a, b], ignore_index=True))):
            f, fp, fn, d = [], [], [], []
            for seed in SEEDS:
                call = X.test_predict(train, test, BASE, seed) >= X.KNOB
                s = X.scores(y, call)
                f.append(s["F1"]); fp.append(s["FP"]); fn.append(s["FN"])
                if label.startswith("A ("):
                    calls_a[seed] = call
                else:
                    d.append(X.bootstrap_diff(test, calls_a[seed], call, n=500, seed=seed))
            dd = np.array(d).mean(axis=0) if d else (0, 0, 0)
            lines.append(f"| {label} | {np.mean(f):.4f} | {np.mean(fp):.1f} | {np.mean(fn):.1f} | "
                         f"{dd[0]:+.4f} ({dd[1]:+.4f}, {dd[2]:+.4f}) |")
            print(lines[-1], flush=True)
        lines.append("")
    with open(os.path.join(HERE, "..", "results", "more_data_exp.md"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
