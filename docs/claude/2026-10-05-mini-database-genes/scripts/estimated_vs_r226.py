#!/usr/bin/env python3
"""protal's conservation factors estimated on a mini world (its database's gene_conservation.tsv) against the real
GTDB r226 factors (gene_rates_r226.tsv): Spearman over all genes and per marker set.

Usage: estimated_vs_r226.py GENE_RATES_R226_TSV DB_FOLDER [DB_FOLDER ...]
"""
import csv
import os
import statistics
import sys


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2
        i = j + 1
    return r


def spearman(a, b):
    return statistics.correlation(ranks(a), ranks(b))


real = {}
with open(sys.argv[1]) as fh:
    for r in csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"):
        real.setdefault(r["marker"], (float(r["within_factor"]), set()))[1].add(r["set"])
for db in sys.argv[2:]:
    ids = dict(line.split() for line in open(os.path.join(db, "gene2geneid.tsv")))
    est = {}
    with open(os.path.join(db, "gene_conservation.tsv")) as fh:
        for line in fh:
            f = line.split("\t")
            if f[0].isdigit():
                est[f[0]] = float(f[1])
    rows = [(m, real[m][0], est[ids[m]], real[m][1]) for m in real if ids.get(m) in est]
    out = [f"{os.path.basename(db)}: {len(rows)} genes, Spearman {spearman([r[1] for r in rows], [r[2] for r in rows]):+.3f}"]
    for s in ("bac120", "ar53"):
        sel = [r for r in rows if r[3] == {s}]
        out.append(f"{s}-only {len(sel)} genes {spearman([r[1] for r in sel], [r[2] for r in sel]):+.3f}")
    print("; ".join(out))
