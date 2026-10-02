#!/usr/bin/env python3
"""How well the r226 v3 training data pin each depth's knob: per half decade of sample depth (log10 of the sample's
fragments over its taxa, as the trainer and protal read it), the threshold with the highest F1 on the scores with
species held out (p_species of trained_model*.predictions.tsv.gz), the range of thresholds within 0.002 F1 of it, and
the 10th-90th percentile of the best threshold over 200 bootstrap resamples of the bin's samples.
Usage: knob_noise.py V3_DIR"""
import os
import sys

import numpy as np
import pandas as pd

GRID = np.round(np.arange(0.05, 0.955, 0.01), 2)
RNG = np.random.default_rng(1)


def f1s(y, p):
    out = []
    for t in GRID:
        call = p >= t
        tp, fp, fn = (call & (y == 1)).sum(), (call & (y == 0)).sum(), (~call & (y == 1)).sum()
        out.append(2 * tp / (2 * tp + fp + fn) if tp + fp + fn else 0.0)
    return np.array(out)


def main():
    v3 = sys.argv[1]
    for t, suffix in (("pe", ""), ("se", "_se"), ("pb", "_pb"), ("ont", "_ont")):
        df = pd.read_csv(os.path.join(v3, f"trained_model{suffix}.predictions.tsv.gz"), sep="\t",
                         usecols=["meta_sample", "truth", "p_species", "fragments"])
        df = df[df["p_species"].notna()]
        df["depth"] = np.log10(np.maximum(df.groupby("meta_sample")["fragments"].transform("sum"), 1.0))
        df["bin"] = np.floor(df["depth"] / 0.5) * 0.5
        print(f"## {t}")
        rows = []
        for b, part in df.groupby("bin"):
            y, p = part["truth"].to_numpy(), part["p_species"].to_numpy()
            f = f1s(y, p)
            best = GRID[int(np.argmax(f))]
            near = GRID[f >= f.max() - 0.002]
            samples = part["meta_sample"].unique()
            by_sample = {s: g for s, g in part.groupby("meta_sample")}
            boot = []
            for _ in range(200):
                pick = pd.concat([by_sample[s] for s in RNG.choice(samples, len(samples))])
                boot.append(GRID[int(np.argmax(f1s(pick["truth"].to_numpy(), pick["p_species"].to_numpy())))])
            rows.append({"bin": f"{b:.1f}-{b + 0.5:.1f}", "samples": len(samples), "taxa": len(part),
                         "present": int(y.sum()), "best knob": best, "F1 there": round(float(f.max()), 4),
                         "F1 at 0.5": round(float(f[GRID == 0.5][0]), 4),
                         "within 0.002": f"{near.min():.2f}-{near.max():.2f}",
                         "bootstrap 10-90%": f"{np.percentile(boot, 10):.2f}-{np.percentile(boot, 90):.2f}"})
        print(pd.DataFrame(rows).to_string(index=False))
        print()


if __name__ == "__main__":
    main()
