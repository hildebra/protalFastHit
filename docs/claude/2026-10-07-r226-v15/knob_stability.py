#!/usr/bin/env python3
"""How stable is the global knob, and what would steadier rules choose? From ablate.py's v15 refits (training rows with
p_species, test rows with the final model's p), per read type:

  - the knob as the trainer chooses it: the threshold of highest weighted F1 with species held out (scenario rows 0.25),
    kept if it gains 0.002 over 0.5;
  - its bootstrap over the training samples (resampled with replacement, 500 times): the best threshold's 5/50/95%
    quantiles, and the share of resamples in which the chosen knob gains over 0.5 at all;
  - rules: "median" (the median of the resamples' best thresholds, kept if it gains over 0.5 in 95% of them), and
    "test" (the median knob, kept only if it does not lose weighted F1 on the test set);
  - each knob's F1 on the design's test rows, the scenarios' hold-out rows, and both weighted as in training.

    python3 knob_stability.py --dir ~/v15/ablate --variant v15 --read-types pe,se,pb,ont
"""
import argparse
import os

import numpy as np
import pandas as pd

GRID = np.round(np.arange(0.30, 0.901, 0.01), 2)


def counts(samples, y, p, w):
    """Per sample (sorted) and threshold: weighted TP, FP, FN. -> array samples x thresholds x 3."""
    order = pd.Index(np.unique(samples))
    idx = order.get_indexer(samples)
    out = np.zeros((len(order), len(GRID), 3))
    for j, t in enumerate(GRID):
        c = p >= t
        out[:, j, 0] = np.bincount(idx, weights=w * (c & y), minlength=len(order))
        out[:, j, 1] = np.bincount(idx, weights=w * (c & ~y), minlength=len(order))
        out[:, j, 2] = np.bincount(idx, weights=w * (~c & y), minlength=len(order))
    return out


def f1(m):
    tp, fp, fn = m[..., 0], m[..., 1], m[..., 2]
    return np.where(tp > 0, 2 * tp / np.maximum(2 * tp + fp + fn, 1e-12), 0.0)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--dir", required=True)
    ap.add_argument("--variant", default="v15")
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    ap.add_argument("--reps", type=int, default=500)
    opts = ap.parse_args()
    rng = np.random.default_rng(1)
    half = int(np.flatnonzero(GRID == 0.5)[0])
    rows = []
    for rt in opts.read_types.split(","):
        path = os.path.join(opts.dir, f"{rt}_{opts.variant}.tsv.gz")
        if not os.path.exists(path):
            continue
        d = pd.read_csv(path, sep="\t", usecols=["part", "meta_sample", "meta_scenario", "truth", "p_species", "p"],
                        dtype={"meta_scenario": str}, keep_default_na=False, na_values=[""])
        scenario = d["meta_scenario"].fillna("design").to_numpy() != "design"
        w = np.where(scenario, 0.25, 1.0)
        y = d["truth"].to_numpy() == 1
        tr = d["part"].to_numpy() == "training"
        c = counts(d["meta_sample"].to_numpy()[tr], y[tr], d["p_species"].to_numpy()[tr], w[tr])
        whole = f1(c.sum(0))
        best = float(GRID[whole.argmax()])
        trainer = best if whole.max() >= whole[half] + 0.002 else 0.5
        n = c.shape[0]
        draws = rng.multinomial(n, np.full(n, 1 / n), size=opts.reps).astype(float)
        boot = f1(np.einsum("bs,stk->btk", draws, c))  # reps x thresholds
        bests = GRID[boot.argmax(1)]
        median = float(np.round(np.median(bests), 2))
        j = int(np.flatnonzero(GRID == median)[0])
        steady = median if (boot[:, j] > boot[:, half]).mean() >= 0.95 else 0.5
        te = ~tr
        test = counts(d["meta_sample"].to_numpy()[te], y[te], d["p"].to_numpy()[te], w[te])
        test_design = counts(d["meta_sample"].to_numpy()[te & ~scenario], y[te & ~scenario], d["p"].to_numpy()[te & ~scenario],
                             np.ones((te & ~scenario).sum()))
        test_scen = counts(d["meta_sample"].to_numpy()[te & scenario], y[te & scenario], d["p"].to_numpy()[te & scenario],
                           np.ones((te & scenario).sum()))
        tw, td, ts = f1(test.sum(0)), f1(test_design.sum(0)), f1(test_scen.sum(0))
        checked = steady if steady == 0.5 or tw[int(np.flatnonzero(GRID == steady)[0])] >= tw[half] else 0.5
        for rule, k in (("0.5", 0.5), ("trainer (0.002 rule)", trainer), ("median, 95% of resamples", steady),
                        ("median, test set not worse", checked)):
            i = int(np.flatnonzero(GRID == k)[0])
            rows.append({"read type": rt, "rule": rule, "knob": k, "train F1 (weighted)": whole[i],
                         "test design F1": td[i], "test scenarios F1": ts[i], "test weighted F1": tw[i]})
        print(f"{rt}: best {best:.2f}; bootstrap best 5/50/95% {np.quantile(bests, 0.05):.2f} / {median:.2f} / "
              f"{np.quantile(bests, 0.95):.2f}; the trainer's knob {trainer:.2f} gains over 0.5 in "
              f"{(boot[:, int(np.flatnonzero(GRID == trainer)[0])] > boot[:, half]).mean():.0%} of resamples; "
              f"test design best {GRID[td.argmax()]:.2f}, scenarios best {GRID[ts.argmax()]:.2f}, weighted best "
              f"{GRID[tw.argmax()]:.2f}", flush=True)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
