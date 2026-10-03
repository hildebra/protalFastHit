#!/usr/bin/env python3
"""The risk of the relatives features: a present species beside a more abundant present congener. How often the
simulated samples have one, and the miss rate of such species with and without the features, by how many times
more fragments the congener has. Usage: congeners.py OUT_DIR TABLES_DIR [rt]"""
import sys

import numpy as np
import pandas as pd

pd.set_option("display.width", 220)
out, tables = sys.argv[1:3]
rt = sys.argv[3] if len(sys.argv) > 3 else ""
TAGS = sys.argv[4].split(",") if len(sys.argv) > 4 else ["base", "rel"]
for name, table, pred, pc in (("species held out", f"training{rt}", "predictions", "p_species"),
                              ("test set", f"test{rt}", "test_predictions", "p")):
    t = pd.read_csv(f"{tables}/{table}.tsv", sep="\t", low_memory=False,
                    usecols=["meta_sample", "taxon", "truth", "fragments"])
    t["truth"] = t["truth"].astype(str).str.lower().isin(("1", "true"))
    calls = {}
    for tag in TAGS:
        p = pd.read_csv(f"{out}/{tag}{rt}.{pred}.tsv.gz", sep="\t", low_memory=False,
                        usecols=["meta_sample", "taxon", pc, "meta_lineage_genus"])
        calls[tag] = p.rename(columns={pc: "p_" + tag})
    df = t.merge(calls[TAGS[0]], on=["meta_sample", "taxon"])
    for tag in TAGS[1:]:
        df = df.merge(calls[tag].drop(columns="meta_lineage_genus"), on=["meta_sample", "taxon"])
    # the most fragments of another PRESENT species of the genus in the sample
    pres = df[df["truth"]]
    key = pres["meta_sample"].astype(str) + "|" + pres["meta_lineage_genus"].astype(str)
    srt = pres.assign(k=key).sort_values("fragments", ascending=False)
    first = srt.groupby("k")["fragments"].transform("first")
    second = srt.groupby("k")["fragments"].transform(lambda s: s.iloc[1] if len(s) > 1 else 0)
    top_other = np.where(srt.groupby("k").cumcount() == 0, second, first)
    pres = pres.assign(top_other=pd.Series(top_other, index=srt.index).reindex(pres.index))
    ratio = (pres["top_other"] + 1) / (pres["fragments"] + 1)
    pres = pres.assign(bin=pd.cut(ratio, [0, 1, 3, 10, 30, 100, 1e12],
                                  labels=["no larger congener", "1-3x", "3-10x", "10-30x", "30-100x", ">100x"],
                                  right=False))
    print(f"\n{name}{rt}: present species by the fragments of the largest other present congener over their own")
    rows = []
    for b, g in pres.groupby("bin", observed=False):
        rows.append({"congener": b, "present": len(g), **{f"missed {tag}": int((g["p_" + tag] < 0.5).sum()) for tag in TAGS}})
    print(pd.DataFrame(rows).to_string(index=False))
