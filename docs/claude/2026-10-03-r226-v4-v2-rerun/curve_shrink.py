"""Test-set F1 of a run's models at knob 0.5, at their knob curve, and at the curve with points whose F1 gain over 0.5
(on their training window, species held out) is below DELTA set back to 0.5; and errors by test depth.
Usage: curve_shrink.py TRAINER_DIR RUN_DIR [TEST_DIR] (TEST_DIR: the test tables, default RUN_DIR/test)"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, sys.argv[1])
import random_forest_cmdline as rf

run = sys.argv[2]
test_dir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(run, "test")


def f1(y, call):
    tp, fp, fn = (call & y).sum(), (call & ~y).sum(), (~call & y).sum()
    return 2 * tp / (2 * tp + fp + fn), int(fp), int(fn)


for t, suffix in (("pe", ""), ("se", "_se"), ("pb", "_pb"), ("ont", "_ont")):
    path = os.path.join(run, f"trained_model{suffix}.metrics.json")
    if not os.path.exists(path):
        continue
    m = json.load(open(path))
    points = [pt for pt in m["depth_knobs"]["points"] if pt["knob"] is not None]
    print(f"## {t} points (log10: knob, gain over 0.5): " + ", ".join(f"{pt['log10 fragments']}: {pt['knob']} {pt['F1 at the knob'] - pt['F1 at 0.5']:+.4f}" for pt in points))
    pred = pd.read_csv(os.path.join(run, f"trained_model{suffix}.test_predictions.tsv.gz"), sep="\t",
                       usecols=["meta_sample", "meta_read_pairs", "taxon", "truth", "p"])
    table = pd.read_csv(os.path.join(test_dir, f"training_data{suffix}.tsv"), sep="\t", usecols=["meta_sample", "taxon", "fragments"])
    df = pred.merge(table, on=["meta_sample", "taxon"])
    assert len(df) == len(pred)
    depth = np.log10(np.maximum(df.groupby("meta_sample")["fragments"].transform("sum").to_numpy(float), 1.0))
    y, p = df["truth"].to_numpy() == 1, df["p"].to_numpy()
    out = [f"0.5: {f1(y, p >= 0.5)[0]:.4f}"]
    calls = {"0.5": p >= 0.5}
    for delta in (0.0, 0.002, 0.005, 0.01):
        curve = [(pt["log10 fragments"], pt["knob"] if pt["F1 at the knob"] - pt["F1 at 0.5"] >= delta else 0.5) for pt in points]
        call = p >= rf.knob_at(curve, depth)
        calls[f"curve {delta}"] = call
        F, fp, fn = f1(y, call)
        out.append(f"curve, gain >= {delta}: {F:.4f} ({fp}+{fn}, {sum(k != 0.5 for _, k in curve)} of {len(curve)} points off 0.5)")
    print(f"## {t}: " + " | ".join(out))
    rows = []
    for d, g in df.assign(i=np.arange(len(df))).groupby("meta_read_pairs"):
        idx = g["i"].to_numpy()
        row = {"depth": d, "samples": g.meta_sample.nunique(), "present": int(y[idx].sum())}
        for name in ("0.5", "curve 0.0", "curve 0.005"):
            _, fp, fn = f1(y[idx], calls[name][idx])
            row[name] = f"{fp}+{fn}"
        rows.append(row)
    print(pd.DataFrame(rows).to_string(index=False))
