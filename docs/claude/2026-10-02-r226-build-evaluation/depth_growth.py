#!/usr/bin/env python3
"""Per design depth: absent and present taxa per sample, the sample's total fragments (what 0.7.3's depth
knobs bin by: digits less one), and the false-positive rate at knob 0.5 (species held out), with a power-law
fit of absent taxa per sample against depth. Usage: depth_growth.py MODEL_LOGS_DIR"""
import sys

import numpy as np
import pandas as pd

logs = sys.argv[1]
for t, suffix in (("pe", ""), ("se", "_se"), ("pb", "_pb"), ("ont", "_ont")):
    df = pd.read_csv(f"{logs}/trained_model{suffix}.predictions.tsv.gz", sep="\t", low_memory=False)
    df["present"] = df["truth"].astype(str).str.lower().isin(("1", "true"))
    per = df.groupby(["meta_read_pairs", "meta_sample"]).agg(
        absent=("present", lambda s: int((~s).sum())), present=("present", "sum"), fragments=("fragments", "sum"),
        fp=("p_species", lambda s: 0)).reset_index()
    fp = df[~df["present"] & (df["p_species"] >= 0.5)].groupby(["meta_read_pairs", "meta_sample"]).size()
    per["fp"] = [int(fp.get((d, s), 0)) for d, s in zip(per["meta_read_pairs"], per["meta_sample"])]
    out = per.groupby("meta_read_pairs").agg(samples=("meta_sample", "count"), absent=("absent", "mean"),
                                             present=("present", "mean"), fragments=("fragments", "median"),
                                             fp=("fp", "mean")).reset_index()
    out["knob_bin"] = np.floor(np.log10(out["fragments"].clip(lower=1))).astype(int).clip(2, 6)
    out["fp_rate_%"] = 100 * out["fp"] / out["absent"]
    slope, intercept = np.polyfit(np.log10(out["meta_read_pairs"]), np.log10(out["absent"]), 1)
    print(f"\n{t}: absent taxa per sample ~ depth^{slope:.2f}")
    print(out.round(2).to_string(index=False))
    if t in ("pe", "se"):
        for depth in (2e6, 1e7, 2e7):
            absent = 10 ** (intercept + slope * np.log10(depth))
            print(f"  extrapolated to {depth:.0e} read pairs: {absent:.0f} absent taxa per sample; at the deepest "
                  f"point's FP rate ({out['fp_rate_%'].iloc[-1]:.1f}%) {absent * out['fp_rate_%'].iloc[-1] / 100:.0f} FP per sample")
