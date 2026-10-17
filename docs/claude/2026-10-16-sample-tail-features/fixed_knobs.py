#!/usr/bin/env python3
"""F1 of both retrained models (as built, with the tail features) at the same fixed knobs, so that a knob the trainer chose
otherwise for one of them does not mix into the features' effect: per read type and set, F1 (TP/FP/FN summed) at each
knob for each model, from their calls files' scores (p).

Usage: fixed_knobs.py OUT_DIR [KNOB ...]   (default knobs 0.5 0.69)
"""
import os
import sys

import pandas as pd

out = sys.argv[1]
knobs = [float(k) for k in sys.argv[2:]] or [0.5, 0.69]
for rt in ("pe", "se"):
    frames = {}
    for name in ("as_built", "with_tail"):
        df = pd.read_csv(os.path.join(out, f"{name}_{rt}.calls.tsv.gz"), sep="\t", keep_default_na=False,
                         usecols=["set", "meta_scenario", "meta_sample", "taxon", "truth", "p"], dtype={"meta_scenario": str})
        frames[name] = df
    print(f"## {rt}")
    print("\t".join(["set"] + [f"{n} @{k:g}" for k in knobs for n in ("as built", "with tail")] + [f"diff @{k:g}" for k in knobs]))
    base = frames["as_built"]
    sets = [("held out, design", (base["set"] == "training") & (base["meta_scenario"] == "")),
            ("held out, scenarios", (base["set"] == "training") & (base["meta_scenario"] != "")),
            ("test, design", (base["set"] == "test") & (base["meta_scenario"] == "")),
            ("test, scenarios", (base["set"] == "test") & (base["meta_scenario"] != ""))]
    for s in ("gut", "moderate", "soil", "soil_shallow"):
        m = (base["set"] == "test") & (base["meta_scenario"] == s)
        if m.any():
            sets.append((f"test, {s}", m))
    for label, mask in sets:
        cells, diffs = [], []
        for k in knobs:
            f1s = []
            for name in ("as_built", "with_tail"):
                d = frames[name][mask.to_numpy()]
                y, c = d["truth"].to_numpy() == 1, d["p"].to_numpy(dtype=float) >= k
                tp, fp, fn = (c & y).sum(), (c & ~y).sum(), (~c & y).sum()
                f1s.append(2 * tp / (2 * tp + fp + fn))
                cells.append(f"{f1s[-1]:.4f} ({tp}/{fp}/{fn})")
            diffs.append(f"{f1s[1] - f1s[0]:+.4f}")
        print("\t".join([label] + cells + diffs))
    print()
