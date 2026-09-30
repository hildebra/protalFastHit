#!/usr/bin/env python3
"""Per species x sample: median depth vs number of multi-allelic genes (MRate2 > 0) in the meta."""
import sys, glob, os, csv, statistics
from collections import defaultdict

for meta in sorted(glob.glob(os.path.join(sys.argv[1], "*.meta.tsv"))):
    sp = os.path.basename(meta)[:-9]
    d = defaultdict(list)
    for r in csv.DictReader(open(meta), delimiter="\t"):
        d[r["sample"]].append(r)
    if not d:
        continue
    print(sp)
    for s in sorted(d):
        rows = d[s]
        dep = statistics.median(float(r["mean_vcov_nonzero"]) for r in rows)
        nbad = sum(1 for r in rows if float(r["multi_rate_vcov2"]) > 0)
        nsites = sum(int(r["multi_allelic"]) for r in rows)
        lowh = sum(1 for r in rows if float(r["hcov"]) < 0.3)
        print(f"  {s:10} genes={len(rows):4} median_depth={dep:7.1f} multi_genes={nbad:4} multi_sites={nsites:5} hcov<0.3: {lowh}")
