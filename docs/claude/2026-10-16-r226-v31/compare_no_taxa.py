#!/usr/bin/env python3
"""Compare the models retrained without sample_log_taxa (retrain_without_taxa.sh) with those retrained as built, on
the same rows: F1 as protal calls (each model at its own knob) per set and scenario with a paired bootstrap over
samples, the moderate/soil/shallow-soil hold-out samples one by one (misses among the species with at most 2 marker
reads), and whether the as-built retraining reproduces the build's own calls.

Usage: compare_no_taxa.py OUT_DIR V31_MODEL_LOGS > no_taxa.txt
  OUT_DIR: retrain_without_taxa.sh's out/ (as_built_pe, no_taxa_pe, as_built_se, no_taxa_se .calls.tsv.gz)
  V31_MODEL_LOGS: the build's model_logs (trained_model.calls.tsv.gz, trained_model_se.calls.tsv.gz)
"""
import argparse
import os

import numpy as np
import pandas as pd

COLUMNS = ["meta_sample", "meta_scenario", "meta_marker_reads", "taxon", "truth", "set", "p", "knob", "call"]
RESAMPLES = 1000
THIN = 2


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
        m = test & (df["meta_scenario"] == s)
        if m.any():
            out.append((f"test, {s} hold-out", m))
        m = train & (df["meta_scenario"] == s)
        if m.any():
            out.append((f"held out, {s} hold-in", m))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out")
    ap.add_argument("model_logs")
    opts = ap.parse_args()
    for rt, suffix in (("pe", ""), ("se", "_se")):
        base = calls(os.path.join(opts.out, f"as_built_{rt}.calls.tsv.gz"))
        alt = calls(os.path.join(opts.out, f"no_taxa_{rt}.calls.tsv.gz"))
        built = calls(os.path.join(opts.model_logs, f"trained_model{suffix}.calls.tsv.gz"))
        key = ["set", "meta_sample", "taxon"]
        m = base.merge(alt[key + ["p", "knob", "call"]], on=key, suffixes=("_b", "_n"), validate="one_to_one")
        assert len(m) == len(base) == len(alt), (rt, len(m), len(base), len(alt))
        r = base.merge(built[key + ["p", "call"]], on=key, suffixes=("", "_built"), validate="one_to_one")
        print(f"## {rt}")
        print(f"as-built retraining against the build's calls: {len(r)} of {len(built)} rows matched; max |p - p_built| "
              f"{(r['p'] - r['p_built']).abs().max():.3g}; calls differ at {int((r['call'] != r['call_built']).sum())} rows")
        print(f"knobs: as built {m['knob_b'].iloc[0]:g}, without sample_log_taxa {m['knob_n'].iloc[0]:g}")
        print()
        print("\t".join(["set", "samples", "present", "as built TP/FP/FN", "no taxa TP/FP/FN", "F1 as built",
                         "F1 no taxa", "no taxa - as built", "95% interval"]))
        for label, mask in groups(m):
            d = m[mask]
            sa, a = per_sample(d, d["call_b"].to_numpy() == 1)
            sb, b = per_sample(d, d["call_n"].to_numpy() == 1)
            fa, fb, diff, lo, hi = paired(a, b)
            ta, tb = a.sum(0).astype(int), b.sum(0).astype(int)
            print("\t".join([label, str(len(sa)), str(int(ta[0] + ta[2])), "/".join(map(str, ta)),
                             "/".join(map(str, tb)), f"{fa:.4f}", f"{fb:.4f}", f"{diff:+.4f}",
                             f"[{lo:+.4f}, {hi:+.4f}]"]))
        print()
        print("hold-out samples of moderate, soil and soil_shallow: misses (FN), and among the present species with "
              f"at most {THIN} marker reads (thin), as built -> without sample_log_taxa; FP likewise")
        print("\t".join(["sample", "present", "FN", "thin", "thin FN", "p thin (median)", "FP"]))
        t = m[(m["set"] == "test") & m["meta_scenario"].isin(["moderate", "soil", "soil_shallow"])]
        for sample, g in t.groupby("meta_sample"):
            p = g[g["truth"] == 1]
            thin = p["marker_reads"] <= THIN
            a = g[g["truth"] == 0]
            print("\t".join([sample, str(len(p)),
                             f"{int((p['call_b'] == 0).sum())} -> {int((p['call_n'] == 0).sum())}",
                             str(int(thin.sum())),
                             f"{int(((p['call_b'] == 0) & thin).sum())} -> {int(((p['call_n'] == 0) & thin).sum())}",
                             f"{p.loc[thin, 'p_b'].median():.3f} -> {p.loc[thin, 'p_n'].median():.3f}"
                             if thin.any() else "-",
                             f"{int(a['call_b'].sum())} -> {int(a['call_n'].sum())}"]))
        print()


if __name__ == "__main__":
    main()
