#!/usr/bin/env python3
"""Every model of the clean evaluation (eval2.sh: the cross-scoring, the feature groups and the depth extrapolation) on
its test set(s): F1, FP and FN at knob 0.5, at the model's knob curve and at its target share of false calls, the misses
of present species beside a present congener of 10 times their fragments or more at each, AP, the F1 at the test set's
best threshold, the cross-validated F1 (species held out); and the extrapolation's errors by depth.
Usage: eval2_summary.py EVAL_DIR OUT_DIR SCRIPTS_DIR"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

E, OUT, SCRIPTS = sys.argv[1:4]
sys.path.insert(0, SCRIPTS)
import random_forest_cmdline as rf  # noqa: E402

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
pd.set_option("display.max_rows", 200)


def minor_congener_mask(df):
    present = df["truth"].to_numpy() == 1
    genus = df["meta_lineage_genus"].astype(str).to_numpy()
    sample = df["meta_sample"].astype(str).to_numpy()
    fragments = df["fragments"].to_numpy(dtype=float)
    values = pd.Series(np.where(present, fragments, -1.0))
    groups = values.groupby([sample, genus])
    top = groups.transform(lambda s: s.nlargest(1).iloc[0]).to_numpy()
    second = groups.transform(lambda s: (s.nlargest(2).tolist() + [-1.0])[1]).to_numpy()
    other = np.where(fragments >= top, second, top)
    return present & (other + 1 >= 10 * (fragments + 1))


tables = {}
rows = []
for path in sorted(glob.glob(os.path.join(E, "*.metrics.json"))):
    name = os.path.basename(path)[:-len(".metrics.json")]
    rt = "_se" if name.endswith("_se") else ""
    stem = name[:-3] if rt else name
    trained, rest = stem.split("_", 1)
    features, test = rest.rsplit("_", 1)
    with open(path) as fh:
        m = json.load(fh)
    key = (test, rt)
    if key not in tables:
        tables[key] = rf.load_table(os.path.join(OUT, test, "test", f"training_data{rt}.tsv"), None)
    pred = pd.read_csv(os.path.join(E, name + ".test_predictions.tsv.gz"), sep="\t", low_memory=False)
    df = tables[key].merge(pred[["meta_sample", "taxon", "p", "meta_lineage_genus"]], on=["meta_sample", "taxon"])
    y, p = df["truth"].to_numpy(), df["p"].to_numpy()
    cs = rf.call_scores(df, p)
    fc = m["false_calls"]
    calls = {"0.5": cs >= 0.5,
             "curve": rf.depth_knob_calls(cs, rf.sample_depths(df), [tuple(c) for c in m["depth_knobs"]["curve"]], 0.5),
             "fdr": rf.FalseCallSamples(df, y, p, [tuple(c) for c in fc["curve"]], fc["prior"]).calls(fc["fdr"])}
    minor = minor_congener_mask(df)
    row = {"reads": rt.strip("_") or "pe", "trained on": trained, "features": features.replace("normalized-adjacency", "na")
           .replace("na-relatives", "na+relatives").replace("na-", "na+"), "tested on": test,
           "CV F1": round(m["evaluation"]["species"]["F1"], 4), "AP": round(m["test"]["this one"]["AP"], 4),
           "best F1": round(m["test_best_threshold"]["F1"], 4), "fdr": fc["fdr"], "minor": int(minor.sum())}
    for mode, c in calls.items():
        tp, fp, fn = int((c & (y == 1)).sum()), int((c & (y == 0)).sum()), int((~c & (y == 1)).sum())
        row[f"F1 {mode}"] = round(2 * tp / max(1, 2 * tp + fp + fn), 4)
        row[f"FP {mode}"], row[f"FN {mode}"] = fp, fn
        row[f"minor missed {mode}"] = int((minor & ~c).sum())
    rows.append(row)
frame = pd.DataFrame(rows).sort_values(["reads", "tested on", "trained on", "features"])
frame.to_csv(os.path.join(E, "summary.tsv"), sep="\t", index=False)
print("=== F1 (and CV F1, AP, best-threshold F1) ===")
print(frame[["reads", "trained on", "features", "tested on", "CV F1", "F1 0.5", "F1 curve", "F1 fdr", "AP", "best F1",
             "fdr"]].to_string(index=False))
for mode in ("0.5", "curve", "fdr"):
    print(f"\n=== errors and minor congeners missed at {mode} ===")
    print(frame[["reads", "trained on", "features", "tested on", f"FP {mode}", f"FN {mode}", "minor",
                 f"minor missed {mode}"]].to_string(index=False))
print("\n=== depth extrapolation (trained on samples of up to 200,000 read pairs) ===")
for path in sorted(glob.glob(os.path.join(E, "shallow_*.metrics.json"))):
    with open(path) as fh:
        m = json.load(fh)
    print(os.path.basename(path)[:-len(".metrics.json")])
    for r in m["test_by_depth"]:
        print(f"  {r['depth']:>8}: present {r['present']:4d} absent {r['absent']:5d} | FN/FP 0.5 {r['FN']:3d}/{r['FP']:3d}, "
              f"curve {r['FN curve']:3d}/{r['FP curve']:3d}, fdr {r['FN fdr']:3d}/{r['FP fdr']:3d}")
