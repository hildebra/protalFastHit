#!/usr/bin/env python3
"""Tables of trace_marker_reads.py's OUT.reads.tsv.gz files: the database species' marker reads (held-out species left
out) by fate, by how much of a gene they hold, by their mean base quality and length, and those of the species with no
row.

Usage: aggregate_marker_reads.py READS.tsv.gz ... > tables.txt
"""
import sys

import pandas as pd

FATES = ["own", "own_lowmapq", "elsewhere", "unaligned", "none"]


def table(df, by, label):
    t = pd.crosstab(df[by], df["fate"]).reindex(columns=FATES, fill_value=0)
    t["reads"] = t.sum(axis=1)
    for f in FATES:
        t[f + " %"] = (100 * t[f] / t["reads"]).round(1)
    print(f"### {label}")
    print(t[["reads"] + [f + " %" for f in FATES]].to_string())
    print()


def main():
    df = pd.concat([pd.read_csv(p, sep="\t", dtype={"taxid": str}) for p in sys.argv[1:]], ignore_index=True)
    db = df[df["heldout"] == 0].copy()
    print(f"{len(df)} marker reads in {df['sample'].nunique()} samples; {len(db)} of database species "
          f"({db['taxid'].nunique()} species-sample pairs: "
          f"{db.groupby(['sample', 'taxid']).ngroups})")
    print()
    db["all"] = "all"
    table(db, "all", "database species' marker reads")
    bins = [0, 100, 200, 500, 1000, 10 ** 9]
    labels = ["50-99", "100-199", "200-499", "500-999", ">=1000"]
    db["overlap bin"] = pd.cut(db["overlap"], bins=bins, labels=labels, right=False)
    db["holds"] = db.apply(lambda r: "whole gene" if r["whole"] else f"part: {r['overlap bin']}", axis=1)
    table(db, "holds", "by the gene they hold (the most bases shared with one marker gene)")
    db["meanq bin"] = pd.cut(db["meanq"], bins=[0, 10, 12, 14, 16, 18, 20, 25, 100], right=False)
    table(db, "meanq bin", "by the read's mean base quality")
    db["length bin"] = pd.cut(db["length"], bins=[0, 500, 1000, 2000, 5000, 10000, 20000, 10 ** 9], right=False)
    table(db, "length bin", "by the read's length")
    whole = db[db["whole"] == 1]
    table(whole, "meanq bin", "reads holding a whole gene, by mean base quality")
    if "error_rate" in db.columns and (db["error_rate"] >= 0).any():
        db["accuracy bin"] = pd.cut(1 - db["error_rate"], bins=[0, 0.80, 0.85, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 1.01],
                                    right=False)
        table(db, "accuracy bin", "by the read's accuracy (1 - the simulator's errors / template length)")
        whole = db[db["whole"] == 1]
        table(whole, "accuracy bin", "reads holding a whole gene, by accuracy")

    # species with no row
    sp = db.groupby(["sample", "taxid"]).agg(row=("row", "max"), reads=("read", "size")).reset_index()
    rowless = sp[sp["row"] == 0]
    print(f"### species with marker reads: {len(sp)} (sample, species) pairs; {len(rowless)} with no row "
          f"({100 * len(rowless) / len(sp):.1f}%)")
    print("their marker reads per species: " + ", ".join(f"{k}: {v}" for k, v in
                                                         rowless["reads"].value_counts().sort_index().items()))
    rl = db.merge(rowless[["sample", "taxid"]], on=["sample", "taxid"])
    table(rl, "holds", "the rowless species' marker reads, by the gene they hold")
    table(rl, "meanq bin", "the rowless species' marker reads, by mean base quality")
    # the species with one marker read: how often is that read lost (P of no row)
    one = sp[sp["reads"] == 1]
    print(f"species with one marker read: {len(one)}, of them rowless {int((one['row'] == 0).sum())} "
          f"({100 * (one['row'] == 0).mean():.1f}%)")
    two = sp[sp["reads"] == 2]
    print(f"species with two marker reads: {len(two)}, of them rowless {int((two['row'] == 0).sum())}")


if __name__ == "__main__":
    main()
