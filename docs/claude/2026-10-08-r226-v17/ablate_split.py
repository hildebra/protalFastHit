#!/usr/bin/env python3
"""The ablation's test rows split: the design's test samples and the scenarios' hold-out samples, F1 at each
variant's knob and at 0.5, from ablate_v17.py's per-row predictions.

    python3 ablate_split.py --out ~/v17/ablate --read-type pe > ablate_split_pe.txt
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd


def f1(y, p, t):
    call = p >= t
    tp, fp, fn = (call & (y == 1)).sum(), (call & (y == 0)).sum(), (~call & (y == 1)).sum()
    return 2 * tp / max(1, 2 * tp + fp + fn), int(fp), int(fn)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--read-type", default="pe")
    opts = ap.parse_args(argv)
    knobs = pd.read_csv(os.path.join(opts.out, "results.tsv"), sep="\t")
    rows = []
    for path in sorted(glob.glob(os.path.join(opts.out, f"{opts.read_type}_*.tsv.gz"))):
        variant = os.path.basename(path)[len(opts.read_type) + 1:-len(".tsv.gz")]
        t = pd.read_csv(path, sep="\t", usecols=["part", "meta_scenario", "truth", "p"])
        t = t[t["part"] == "test"]
        knob = float(knobs.loc[(knobs["variant"] == variant) & (knobs["read type"] == opts.read_type), "knob"].iloc[0])
        y = t["truth"].to_numpy()
        sc = t["meta_scenario"].fillna("").astype(str).to_numpy()
        for label, sel in (("design test", sc == ""), ("scenario hold-out", sc != ""),
                           *(((s, sc == s) for s in sorted(set(sc) - {""})))):
            if not sel.any():
                continue
            a, fp, fn = f1(y[sel], t["p"].to_numpy()[sel], knob)
            b, _, _ = f1(y[sel], t["p"].to_numpy()[sel], 0.5)
            rows.append({"variant": variant, "rows": label, "n": int(sel.sum()), "knob": knob, "F1 at knob": round(a, 4),
                         "FP/FN": f"{fp}/{fn}", "F1 at 0.5": round(b, 4)})
    out = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(out.pivot(index="variant", columns="rows", values="F1 at knob").to_string())
    print()
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
