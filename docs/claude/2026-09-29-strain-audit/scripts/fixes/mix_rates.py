#!/usr/bin/env python3
"""Pooled multi-allelic rate per sample (sum multi_allelic / sum counts_vcov2) of a .meta.tsv.
usage: mix_rates.py META [META...]"""
import sys
from collections import defaultdict

for path in sys.argv[1:]:
    multi, pos, genes, bad = defaultdict(int), defaultdict(int), defaultdict(int), defaultdict(int)
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            s = r["sample"]
            multi[s] += int(r["multi_allelic"])
            pos[s] += int(r["counts_vcov2"])
            genes[s] += 1
            bad[s] += int(r["multi_allelic"]) > 0
    print(f"== {path}")
    for s in sorted(multi, key=lambda s: -multi[s] / max(pos[s], 1)):
        print(f"  {s:12s} rate {multi[s] / max(pos[s], 1):.5f}  multi {multi[s]:5d}  positions {pos[s]:7d}  multi-allelic genes {bad[s]}/{genes[s]}")
