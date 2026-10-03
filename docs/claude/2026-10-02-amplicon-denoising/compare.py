#!/usr/bin/env python3
"""base (normalized+adjacency) against rel (+ genus_skew, family_skew, genus_share, genus_spill) on the r226 v3
tables: species held out and the independent test set, at knob 0.5 and at each model's own knob curve, by depth
and by the source of the false positives. Usage: compare.py OUT_DIR TABLES_DIR"""
import re
import sys

import numpy as np
import pandas as pd

pd.set_option("display.width", 220)
pd.set_option("display.max_columns", 30)
out, tables = sys.argv[1:3]
TAGS = sys.argv[3].split(",") if len(sys.argv) > 3 else ["base", "rel"]


def curve_of(report):
    m = re.search(r"F1 of species held out: .* at the knob curve \(([^)]*)\)", open(report).read())
    return [tuple(map(float, p.split(":"))) for p in m.group(1).split(",")] if m else None


def f1(y, c):
    tp, fp, fn = int((y & c).sum()), int((~y & c).sum()), int((y & ~c).sum())
    return 2 * tp / max(1, 2 * tp + fp + fn), fp, fn


def calls(df, p, curve):
    depth = np.log10(np.maximum(df.groupby("meta_sample")["fragments"].transform("sum").to_numpy(), 1))
    xs, ks = zip(*curve)
    return p >= np.interp(depth, xs, ks)


rows, depth_rows, source_rows = [], [], []
for rt in ("", "_se", "_pb", "_ont"):
    for tag in TAGS:
        curve = curve_of(f"{out}/{tag}{rt}.report.txt")
        sp = pd.read_csv(f"{out}/{tag}{rt}.predictions.tsv.gz", sep="\t", low_memory=False)
        te = pd.read_csv(f"{out}/{tag}{rt}.test_predictions.tsv.gz", sep="\t", low_memory=False)
        tt = pd.read_csv(f"{tables}/test{rt}.tsv", sep="\t", low_memory=False, usecols=["meta_sample", "taxon", "fragments"])
        te = te.merge(tt, on=["meta_sample", "taxon"], how="left")
        for name, df, pc in (("species held out", sp, "p_species"), ("test set", te, "p")):
            y = df["truth"].astype(str).str.lower().isin(("1", "true")).to_numpy()
            p = df[pc].to_numpy()
            c05, cc = p >= 0.5, calls(df, p, curve)
            a, b = f1(y, c05), f1(y, cc)
            rows.append({"read type": rt.strip("_") or "pe", "model": tag, "evaluated on": name,
                         "F1 0.5": round(a[0], 4), "FP 0.5": a[1], "FN 0.5": a[2],
                         "F1 curve": round(b[0], 4), "FP curve": b[1], "FN curve": b[2]})
            novel = df["meta_novel_level"].fillna("").astype(str).replace("nan", "")
            src = np.where(novel != "", "near held-out", "near db species")
            for s in ("near held-out", "near db species"):
                m = src == s
                source_rows.append({"read type": rt.strip("_") or "pe", "model": tag, "evaluated on": name,
                                    "source": s, "FP 0.5": int((~y & c05 & m).sum()),
                                    "FP curve": int((~y & cc & m).sum())})
            for d, g in df.assign(_y=y, _c05=c05, _cc=cc).groupby("meta_read_pairs"):
                depth_rows.append({"read type": rt.strip("_") or "pe", "evaluated on": name, "depth": d, "model": tag,
                                   "FP 0.5": int((~g._y & g._c05).sum()), "FN 0.5": int((g._y & ~g._c05).sum()),
                                   "FP curve": int((~g._y & g._cc).sum()), "FN curve": int((g._y & ~g._cc).sum())})
        print(f"{tag}{rt} knob curve: " + ",".join(f"{x:.2f}:{k:.2f}" for x, k in curve))

print("\n" + pd.DataFrame(rows).to_string(index=False))
print("\nFalse positives by source (near held-out: a species or clade the database lacks is nearest)")
print(pd.DataFrame(source_rows).pivot_table(index=["read type", "evaluated on", "source"], columns="model",
                                            values=["FP 0.5", "FP curve"]).to_string())
print("\nBy depth (paired-end)")
d = pd.DataFrame(depth_rows)
print(d[d["read type"] == "pe"].pivot_table(index=["evaluated on", "depth"], columns="model",
                                             values=["FP 0.5", "FN 0.5", "FP curve", "FN curve"]).to_string())
