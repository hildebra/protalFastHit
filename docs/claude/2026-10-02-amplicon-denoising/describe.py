#!/usr/bin/env python3
"""The r226 v3 paired-end calls of the base model (at knob 0.5: species held out, or the test set) by the
taxon's skew to the most abundant relative in its sample: are the false positives the low-skew neighbours of an
abundant relative (the UNOISE/DADA2 case), and how many real species sit in the same place?
Usage: describe.py TABLE PREDICTIONS [P_COLUMN]
genus_skew = log10((fragments + 1) / (most fragments of a congener in the sample + 1)): below 0 a congener has more;
above 0 the taxon is the largest of its genus there, or alone."""
import sys

import numpy as np
import pandas as pd

pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 30)
table, preds = sys.argv[1:3]
PC = sys.argv[3] if len(sys.argv) > 3 else "p_species"
t = pd.read_csv(table, sep="\t", low_memory=False)
p = pd.read_csv(preds, sep="\t", low_memory=False, usecols=["meta_sample", "taxon", PC]).rename(columns={PC: "p_species"})
t["truth"] = t["truth"].astype(str).str.lower().isin(("1", "true"))
df = t.merge(p, on=["meta_sample", "taxon"], how="inner")
print(f"{len(df)} of {len(t)} rows joined")
df["call"] = df["p_species"] >= 0.5
novel = df["meta_novel_level"].fillna("").astype(str).replace("nan", "")
df["source"] = np.where(novel != "", "near held-out", "near db species")
df["skew_bin"] = pd.cut(df["genus_skew"], [-9, -3, -2, -1, 0, 9],
                        labels=["<-3", "-3..-2", "-2..-1", "-1..0 (incl. ties)", ">0 (largest or alone)"])
fp = df[df["call"] & ~df["truth"]]
tp = df[df["call"] & df["truth"]]
fn = df[~df["call"] & df["truth"]]
print(f"\nTP {len(tp)} FP {len(fp)} FN {len(fn)}")
print("\nFP by genus_skew bin x source:")
print(pd.crosstab(fp["skew_bin"], fp["source"], margins=True))
print("\nTP by genus_skew bin; FN by genus_skew bin:")
print(pd.concat([tp["skew_bin"].value_counts().rename("TP"), fn["skew_bin"].value_counts().rename("FN"),
                 fp["skew_bin"].value_counts().rename("FP")], axis=1).sort_index())
print("\nsingle-fragment calls: by whether a congener in the sample has >=10 / >=100 fragments")
one = df[df["call"] & (df["fragments"] <= 1)]
for thr in (10, 100, 1000):
    near = one["genus_top"] >= thr
    print(f"  congener >= {thr:5d}: FP {int((near & ~one['truth']).sum()):4d}  TP {int((near & one['truth']).sum()):4d} | "
          f"no such congener: FP {int((~near & ~one['truth']).sum()):4d}  TP {int((~near & one['truth']).sum()):4d}")
print("\nby depth: FP and TP with genus_skew < -2 (a congener >100x more fragments)")
low = df["genus_skew"] < -2
print(df.assign(low=low).groupby("meta_read_pairs").apply(
    lambda g: pd.Series({"FP": int((g.call & ~g.truth).sum()), "FP_low": int((g.call & ~g.truth & g.low).sum()),
                         "TP": int((g.call & g.truth).sum()), "TP_low": int((g.call & g.truth & g.low).sum()),
                         "absent_low": int((~g.truth & g.low).sum()), "present_low": int((g.truth & g.low).sum())}),
    include_groups=False))
