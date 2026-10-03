#!/usr/bin/env python3
"""Summarize the sample-depth / genus-size experiment: per model, held-out and test F1 at knob 0.5, at the knob curve and
at the best test threshold (from the trainer logs), plus test false positives and misses by fragments bin at 0.5 and at
the curve (from the test predictions). Usage: summarize_depth.py OUTDIR TABLEDIR"""
import glob, os, re, sys
import numpy as np
import pandas as pd
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
out, tables = sys.argv[1], sys.argv[2]

def grab(text, pattern, default=np.nan, cast=float):
    m = re.search(pattern, text)
    return cast(m.group(1)) if m else default

def curve_of(text):
    m = re.search(r"knob curve \(([^)]*)\)", text)
    if not m:
        return None
    return [tuple(float(x) for x in p.split(":")) for p in m.group(1).split(",")]

rows, bins = [], []
for log in sorted(glob.glob(f"{out}/*.log")):
    tag = os.path.basename(log)[:-4]
    rt = "_" + tag.split("_")[-1] if tag.endswith(("_se", "_pb", "_ont")) else ""
    text = open(log).read()
    if "ALLDONE" in text or "## Independent test set" not in text:
        continue
    held = grab(text, r"F1 of species held out: ([0-9.]+) at 0\.5")
    held_curve = grab(text, r"F1 of species held out: [0-9.]+ at 0\.5, ([0-9.]+) at the knob curve")
    test_part = text.split("## Independent test set")[1]
    test05 = grab(test_part, r"at knob 0\.5: F1 ([0-9.]+)")
    test_curve = grab(test_part, r"as protal calls by default: F1 ([0-9.]+)")
    fp_curve = grab(test_part, r"as protal calls by default: F1 [0-9.]+, ([0-9]+) false positives", cast=int)
    fn_curve = grab(test_part, r"false positives, ([0-9]+) false negatives", cast=int)
    fp05 = grab(test_part, r"at knob 0\.5: F1 [0-9.]+, ([0-9]+) and", cast=int)
    fn05 = grab(test_part, r"at knob 0\.5: F1 [0-9.]+, [0-9]+ and ([0-9]+)", cast=int)
    best = grab(test_part, r"highest F1 on the test set at threshold [0-9.]+ \(F1 ([0-9.]+)")
    logloss = grab(test_part, r"this one\s+\d+\s+\d+\s+[0-9.]+\s+[0-9.]+\s+([0-9.]+)")
    ap = grab(test_part, r"this one\s+\d+\s+\d+\s+[0-9.]+\s+([0-9.]+)")
    rows.append(dict(model=tag, held_out_05=held, held_out_curve=held_curve, test_05=test05, test_curve=test_curve, test_best=best,
                     test_AP=ap, test_logloss=logloss, FP_05=fp05, FN_05=fn05, FP_curve=fp_curve, FN_curve=fn_curve))
    # by fragments bin on the test set
    pred = f"{out}/{tag}.test_predictions.tsv.gz"
    if os.path.exists(pred):
        p = pd.read_csv(pred, sep="\t", low_memory=False)
        if "fragments" not in p.columns:
            t = pd.read_csv(f"{tables}/test{rt}.tsv", sep="\t", usecols=["meta_sample", "taxon", "fragments"], low_memory=False)
            p = p.merge(t, on=["meta_sample", "taxon"], how="left")
        curve = curve_of(text)
        depth = np.log10(np.maximum(p.groupby("meta_sample").fragments.transform("sum").astype(float), 1))
        knob = np.interp(depth, *zip(*curve)) if curve else 0.5
        p["call05"] = p.p >= 0.5
        p["callc"] = p.p >= knob
        fb = pd.cut(p.fragments, [0, 1, 2, 9, 99, 1e12], labels=["1", "2", "3-9", "10-99", ">=100"])
        for name, call in (("0.5", p.call05), ("curve", p.callc)):
            g = pd.DataFrame({"bin": fb, "FP": (p.truth == 0) & call, "FN": (p.truth == 1) & ~call}).groupby("bin", observed=False).sum()
            bins.append(dict(model=tag, at=name, **{f"FP_{b}": int(g.FP[b]) for b in g.index}, **{f"FN_{b}": int(g.FN[b]) for b in g.index}))
print(pd.DataFrame(rows).set_index("model").round(4).to_string())
print()
print(pd.DataFrame(bins).set_index(["model", "at"]).to_string())
