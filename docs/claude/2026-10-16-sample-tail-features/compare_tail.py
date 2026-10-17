#!/usr/bin/env python3
"""Compare the models retrained with the sample tail features (retrain_with_tail.sh: the default set + tail) with those
retrained as built, on the same rows: F1 as protal calls (each model at its own knob) per set and scenario with a paired
bootstrap over samples; the moderate/soil/shallow-soil hold-out samples one by one (misses among the species with at
most 2 marker reads, their median score, false positives); and the tail features' importance.

Usage: compare_tail.py OUT_DIR > tail.txt
  OUT_DIR: retrain_with_tail.sh's out/ (as_built_pe, with_tail_pe, as_built_se, with_tail_se .calls.tsv.gz and
  .varimp.tsv)
"""
import argparse
import os

import numpy as np
import pandas as pd

COLUMNS = ["meta_sample", "meta_scenario", "meta_marker_reads", "taxon", "truth", "set", "p", "knob", "call"]
RESAMPLES = 1000
THIN = 2
TAIL = ("sample_thin_share", "sample_thin_close_share")


def calls(path):
    df = pd.read_csv(path, sep="\t", usecols=COLUMNS, dtype={"meta_scenario": str, "meta_sample": str},
                     keep_default_na=False)
    df["p"] = df["p"].astype(float)
    df["marker_reads"] = pd.to_numeric(df["meta_marker_reads"], errors="coerce")
    return df


def per_sample(df, call):
    y = df["truth"].to_numpy() == 1
    t = pd.DataFrame({"s": df["meta_sample"].to_numpy(), "tp": call & y, "fp": call & ~y, "fn": ~call & y})
    g = t.groupby("s", sort=True)[["tp", "fp", "fn"]].sum()
    return g.index.to_numpy(), g.to_numpy(dtype=float, copy=True)


def f1(c):
    tp, fp, fn = c[..., 0], c[..., 1], c[..., 2]
    return 2 * tp / np.maximum(2 * tp + fp + fn, 1)


def paired(a, b, seed=1):
    rng = np.random.default_rng(seed)
    n = len(a)
    w = rng.multinomial(n, np.full(n, 1 / n), size=RESAMPLES).astype(float)
    d = f1(w @ b) - f1(w @ a)
    return f1(a.sum(0)), f1(b.sum(0)), f1(b.sum(0)) - f1(a.sum(0)), np.quantile(d, 0.025), np.quantile(d, 0.975)


def groups(df):
    train, test, scen = df["set"] == "training", df["set"] == "test", df["meta_scenario"] != ""
    out = [("held out, design", train & ~scen), ("held out, scenarios (hold-in)", train & scen),
           ("test, design", test & ~scen), ("test, scenarios (hold-out)", test & scen)]
    for s in ("gut", "moderate", "soil", "soil_shallow", "host"):
        for label, base in ((f"test, {s} hold-out", test), (f"held out, {s} hold-in", train)):
            m = base & (df["meta_scenario"] == s)
            if m.any():
                out.append((label, m))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    opts = ap.parse_args()
    for rt in ("pe", "se"):
        base = calls(os.path.join(opts.out, f"as_built_{rt}.calls.tsv.gz"))
        alt = calls(os.path.join(opts.out, f"with_tail_{rt}.calls.tsv.gz"))
        key = ["set", "meta_sample", "taxon"]
        m = base.merge(alt[key + ["p", "knob", "call"]], on=key, suffixes=("_b", "_t"), validate="one_to_one")
        assert len(m) == len(base) == len(alt), (rt, len(m), len(base), len(alt))
        print(f"## {rt}")
        print(f"knobs: as built {m['knob_b'].iloc[0]:g}, with the tail features {m['knob_t'].iloc[0]:g}")
        varimp = pd.read_csv(os.path.join(opts.out, f"with_tail_{rt}.varimp.tsv"), sep="\t")
        name, value = varimp.columns[0], varimp.columns[-1]
        ranked = varimp.sort_values(value, ascending=False).reset_index(drop=True)
        for feature in TAIL:
            row = ranked[ranked[name] == feature]
            if len(row):
                print(f"{feature}: importance {float(row[value].iloc[0]):.4f}, rank {int(row.index[0]) + 1} of {len(ranked)}")
        print()
        print("\t".join(["set", "samples", "present", "as built TP/FP/FN", "with tail TP/FP/FN", "F1 as built",
                         "F1 with tail", "with tail - as built", "95% interval"]))
        for label, mask in groups(m):
            d = m[mask]
            sa, a = per_sample(d, d["call_b"].to_numpy() == 1)
            sb, b = per_sample(d, d["call_t"].to_numpy() == 1)
            fa, fb, diff, lo, hi = paired(a, b)
            ta, tb = a.sum(0).astype(int), b.sum(0).astype(int)
            print("\t".join([label, str(len(sa)), str(int(ta[0] + ta[2])), "/".join(map(str, ta)),
                             "/".join(map(str, tb)), f"{fa:.4f}", f"{fb:.4f}", f"{diff:+.4f}",
                             f"[{lo:+.4f}, {hi:+.4f}]"]))
        print()
        print(f"hold-out samples of moderate, soil and soil_shallow: misses (FN), and among the present species with at "
              f"most {THIN} marker reads (thin), as built -> with the tail features; FP likewise")
        print("\t".join(["sample", "present", "FN", "thin", "thin FN", "p thin (median)", "FP"]))
        t = m[(m["set"] == "test") & m["meta_scenario"].isin(["moderate", "soil", "soil_shallow"])]
        for sample, g in sorted(t.groupby("meta_sample")):
            p = g[g["truth"] == 1]
            thin = p["marker_reads"] <= THIN
            a = g[g["truth"] == 0]
            print("\t".join([sample, str(len(p)),
                             f"{int((p['call_b'] == 0).sum())} -> {int((p['call_t'] == 0).sum())}",
                             str(int(thin.sum())),
                             f"{int(((p['call_b'] == 0) & thin).sum())} -> {int(((p['call_t'] == 0) & thin).sum())}",
                             f"{p.loc[thin, 'p_b'].median():.3f} -> {p.loc[thin, 'p_t'].median():.3f}"
                             if thin.any() else "-",
                             f"{int(a['call_b'].sum())} -> {int(a['call_t'].sum())}"]))
        print()


if __name__ == "__main__":
    main()
