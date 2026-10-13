#!/usr/bin/env python3
"""Some taxa's features side by side across protal training dumps (<profile>.truth_annotated), with each sample's
median over its true positives.

usage: compare_rows.py FEATURE[,FEATURE...] DUMP:TAXID[,TAXID...] [DUMP:TAXID ...]
"""
import os
import sys

import pandas as pd

features = sys.argv[1].split(",")
columns = []
for spec in sys.argv[2:]:
    dump, taxids = spec.rsplit(":", 1)
    table = pd.read_csv(dump, sep="\t")
    sample = os.path.basename(dump).split(".profile")[0]
    for taxid in taxids.split(","):
        row = table[table["taxon"].astype(str) == taxid].iloc[0]
        columns.append((f"{sample}:{row['taxon_name'].split('__')[-1]}", row))
    tp = table[(table["truth"] == 1) & (table["prediction"] == 1)]
    present = [c for c in features + ["probability"] if c in tp]
    columns.append((f"{sample}:TP median (n={len(tp)})", tp[present].median()))
print("feature\t" + "\t".join(name for name, _ in columns))
for f in ["probability"] + features:
    print(f + "\t" + "\t".join(f"{float(r[f]):.4g}" if f in r.index else "NA" for _, r in columns))
