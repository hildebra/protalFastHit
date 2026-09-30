#!/usr/bin/env python3
"""Which rows and columns differ between two training tables."""
import sys
import pandas as pd

a = pd.read_csv(sys.argv[1], sep="\t", float_precision="round_trip")
b = pd.read_csv(sys.argv[2], sep="\t", float_precision="round_trip")
print("shapes", a.shape, b.shape)
key = ["meta_sample", "taxon"]
m = a.merge(b, on=key, how="outer", suffixes=("_a", "_b"), indicator=True)
print(m["_merge"].value_counts().to_dict())
both = m[m["_merge"] == "both"]
cols = [c for c in a.columns if c not in key]
for c in cols:
    x, y = both[c + "_a"], both[c + "_b"]
    diff = ~((x == y) | (x.isna() & y.isna()))
    if diff.any():
        print(f"{c}: {int(diff.sum())} rows differ, e.g.", list(zip(x[diff].head(3), y[diff].head(3))))
print("samples with differing rows:", sorted(set(both.loc[
    ~((both["fragments_a"] == both["fragments_b"])), "meta_sample"])))
