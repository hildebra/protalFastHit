#!/usr/bin/env python3
"""Two builds' presence models on the same rows, from the builds' own predictions (no refit), with a paired bootstrap
over samples for the difference in F1.

The r226 v14 and v15 builds profiled the same samples (compare_tables.py: same rows, same features but v15's three
complexity columns), so every difference here is the model's. Sets, each at each build's knob (metrics.json
global_knob, else 0.5), as protal calls:

  design test                        the design's test samples, final model (*.test_predictions.tsv.gz, p)
  <scenario> hold-out                the scenario's test samples, final model (*.scenario_predictions.tsv.gz, p),
                                     where both builds have the file
  <scenario> training, species/samples held out
                                     the training rows (*.predictions.tsv.gz): p_species (5 folds by taxon),
                                     p_samples (5 folds by sample); "design" is the design's training samples

The 95% interval of B - A resamples the set's samples with replacement (2000 times), both builds' calls together.

    python3 compare_builds.py --builds v14=local/v14,v15=local/v15 > compare_builds.txt
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

SUFFIX = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
COLS = ["meta_sample", "meta_scenario", "taxon", "truth"]


def knob_of(build, rt):
    with open(os.path.join(build, "model_logs", f"trained_model{SUFFIX[rt]}.metrics.json")) as fh:
        k = (json.load(fh).get("depth_knobs") or {}).get("global_knob")
    return float(k) if k is not None else 0.5


def read(build, rt, kind, score):
    path = os.path.join(build, "model_logs", f"trained_model{SUFFIX[rt]}.{kind}.tsv.gz")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", usecols=COLS + [score], dtype={"meta_scenario": str, "meta_sample": str})
    df["meta_scenario"] = df["meta_scenario"].fillna("design")
    return df.rename(columns={score: "p"})


def counts(df, knob):
    """Per sample: TP, FP, FN at the knob (arrays in the order of the sorted sample names)."""
    y, c = df["truth"].to_numpy() == 1, df["p"].to_numpy() >= knob
    g = pd.DataFrame({"s": df["meta_sample"].to_numpy(), "tp": c & y, "fp": c & ~y, "fn": ~c & y})
    t = g.groupby("s")[["tp", "fp", "fn"]].sum().sort_index()
    return t.index.to_numpy(), t.to_numpy(dtype=float)


def f1(m):
    tp, fp, fn = m[..., 0], m[..., 1], m[..., 2]
    return np.where(tp > 0, 2 * tp / np.maximum(2 * tp + fp + fn, 1), 0.0)


def compare(a, b, ka, kb, rng, reps):
    if len(a) != len(b) or not ((a["meta_sample"].to_numpy() == b["meta_sample"].to_numpy()).all()
                                and (a["taxon"].to_numpy() == b["taxon"].to_numpy()).all()):
        raise SystemExit("the builds' rows differ: compare_builds.py needs the same rows in the same order")
    sa, ca = counts(a, ka)
    sb, cb = counts(b, kb)
    assert (sa == sb).all()
    fa, fb = float(f1(ca.sum(0))), float(f1(cb.sum(0)))
    draws = rng.multinomial(len(sa), np.full(len(sa), 1 / len(sa)), size=reps).astype(float)
    diff = f1(draws @ cb) - f1(draws @ ca)
    lo, hi = np.percentile(diff, [2.5, 97.5])
    return {"rows": len(a), "samples": len(sa), "F1 A": fa, "F1 B": fb, "B - A": fb - fa, "95% low": lo, "95% high": hi,
            "FP A": int(ca[:, 1].sum()), "FP B": int(cb[:, 1].sum()), "FN A": int(ca[:, 2].sum()),
            "FN B": int(cb[:, 2].sum())}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--builds", required=True, help="A=DIR,B=DIR")
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    opts = ap.parse_args()
    (na, da), (nb, db) = [s.split("=", 1) for s in opts.builds.split(",")]
    rng = np.random.default_rng(opts.seed)
    rows = []
    for rt in opts.read_types.split(","):
        ka, kb = knob_of(da, rt), knob_of(db, rt)
        pairs = [("design test", read(da, rt, "test_predictions", "p"), read(db, rt, "test_predictions", "p"))]
        ha, hb = read(da, rt, "scenario_predictions", "p"), read(db, rt, "scenario_predictions", "p")
        if ha is not None and hb is not None:
            for sc in sorted(set(ha["meta_scenario"])):
                pairs.append((f"{sc} hold-out", ha[ha["meta_scenario"] == sc].reset_index(drop=True),
                              hb[hb["meta_scenario"] == sc].reset_index(drop=True)))
        for score, label in (("p_species", "species held out"), ("p_samples", "samples held out")):
            ta, tb = read(da, rt, "predictions", score), read(db, rt, "predictions", score)
            for sc in ["design"] + sorted(set(ta["meta_scenario"]) - {"design"}):
                pairs.append((f"{sc} training, {label}", ta[ta["meta_scenario"] == sc].reset_index(drop=True),
                              tb[tb["meta_scenario"] == sc].reset_index(drop=True)))
        for name, a, b in pairs:
            r = compare(a, b, ka, kb, rng, opts.reps)
            rows.append({"read type": rt, "set": name, f"knob {na}": ka, f"knob {nb}": kb, **r})
    out = pd.DataFrame(rows).rename(columns={"F1 A": f"F1 {na}", "F1 B": f"F1 {nb}", "B - A": f"{nb} - {na}",
                                             "FP A": f"FP {na}", "FP B": f"FP {nb}", "FN A": f"FN {na}",
                                             "FN B": f"FN {nb}"})
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
