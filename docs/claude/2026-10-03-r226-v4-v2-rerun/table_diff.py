"""Which columns differ between two training tables, on the rows (meta_sample, taxon) both have. Usage: table_diff.py A B"""
import sys

import numpy as np
import pandas as pd

a, b = (pd.read_csv(p, sep="\t") for p in sys.argv[1:3])
key = ["meta_sample", "taxon"]
print(f"rows: {len(a)} and {len(b)}; samples {a.meta_sample.nunique()} and {b.meta_sample.nunique()}")
m = a.merge(b, on=key, suffixes=("_a", "_b"))
print(f"rows on both keys: {len(m)}; only in A {len(a) - len(m)}, only in B {len(b) - len(m)}")
only_a = a.merge(b[key], on=key, how="left", indicator=True).query("_merge == 'left_only'")
print("samples with rows only in A:", only_a.meta_sample.value_counts().head(8).to_dict())
for c in a.columns:
    if c in key or c + "_b" not in m:
        continue
    x, y = m[c + "_a"], m[c + "_b"]
    if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
        d = ~np.isclose(x, y, rtol=1e-9, atol=0, equal_nan=True)
        if d.any():
            rel = (np.abs(x - y) / np.maximum(np.abs(x), 1e-300))[d]
            print(f"{c}: {int(d.sum())} rows differ, largest relative {rel.max():.3g}, median {np.median(rel):.3g}")
    else:
        d = x.fillna("").astype(str) != y.fillna("").astype(str)
        if d.any():
            print(f"{c}: {int(d.sum())} rows differ (text), e.g. {x[d].iloc[0]!r} vs {y[d].iloc[0]!r}")
