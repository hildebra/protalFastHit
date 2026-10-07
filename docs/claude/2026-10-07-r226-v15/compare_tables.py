#!/usr/bin/env python3
"""Are two builds' training tables the same rows with the same feature values? Per read type and set (training,
test): rows, keys (meta_sample, taxon, truth) in order, and for every shared column the rows that differ; the columns
only one build has.

    python3 compare_tables.py --builds v14=local/v14,v15=local/v15
"""
import argparse
import os

import numpy as np
import pandas as pd

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
KEYS = ["meta_sample", "taxon", "truth"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--builds", required=True, help="A=DIR,B=DIR")
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    opts = ap.parse_args()
    (na, da), (nb, db) = [s.split("=", 1) for s in opts.builds.split(",")]
    for rt in opts.read_types.split(","):
        for part in ("training", "test"):
            a = pd.read_csv(os.path.join(da, part, TABLES[rt]), sep="\t", low_memory=False)
            b = pd.read_csv(os.path.join(db, part, TABLES[rt]), sep="\t", low_memory=False)
            only_a = [c for c in a.columns if c not in b.columns]
            only_b = [c for c in b.columns if c not in a.columns]
            print(f"## {rt} {part}: {na} {len(a)} rows, {nb} {len(b)} rows; only {na}: {only_a or '-'}; "
                  f"only {nb}: {only_b or '-'}")
            if len(a) != len(b):
                continue
            same_keys = all((a[k].astype(str).to_numpy() == b[k].astype(str).to_numpy()).all() for k in KEYS)
            print(f"   keys (meta_sample, taxon, truth) in the same order: {same_keys}")
            if not same_keys:
                continue
            differ = []
            for c in a.columns:
                if c not in b.columns or c in KEYS:
                    continue
                x, y = a[c], b[c]
                if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
                    xv, yv = x.to_numpy(float), y.to_numpy(float)
                    bad = ~((xv == yv) | (np.isnan(xv) & np.isnan(yv)))
                    if bad.any():
                        differ.append(f"{c}: {int(bad.sum())} rows, max |diff| "
                                      f"{np.nanmax(np.abs(xv[bad] - yv[bad])):.3g}")
                else:
                    # pandas 3 keeps NaN through astype(str), and NaN != NaN
                    bad = x.fillna("").astype(str).to_numpy() != y.fillna("").astype(str).to_numpy()
                    if bad.any():
                        differ.append(f"{c}: {int(bad.sum())} rows (text)")
            print("   shared columns that differ: " + ("; ".join(differ) if differ else "none"))


if __name__ == "__main__":
    main()
