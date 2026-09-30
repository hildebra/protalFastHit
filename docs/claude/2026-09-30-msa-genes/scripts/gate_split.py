#!/usr/bin/env python3
"""Genes qcmsa's coverage gate removed, by the species' number of MSA samples and their depth; and
why species lost their filtered MSA. usage: gate_split.py STRAINS_DIR"""
import csv
import os
import sys
from collections import defaultdict

d = sys.argv[1]
by_n = defaultdict(lambda: [0, 0, 0])  # samples -> [species, genes removed by coverage, genes in raw]
depths = defaultdict(list)
with open(os.path.join(d, "species.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        sp, n = r["species"], int(r["samples"])
        path = os.path.join(d, sp + ".qcmsa_summary.tsv")
        if r["raw_msa"] == "-" or not os.path.exists(path):
            continue
        counts, status = {}, ""
        for line in open(path):
            f = line.rstrip("\n").split("\t")
            if f[0] == "count":
                counts[f[1]] = int(f[2])
            elif f[0] == "status":
                status = f[1] + (" (" + f[3] + ")" if len(f) > 3 and f[3] else "")
        if status.startswith("no_msa"):
            print(f"{sp}: {status}")
        k = n if n < 5 else 5
        by_n[k][0] += 1
        by_n[k][1] += counts.get("genes_filtered_coverage", 0)
        by_n[k][2] += counts.get("genes_in", 0)
        # mean depth per sample of the species (vertical_coverage over its genes)
        per = defaultdict(list)
        with open(os.path.join(d, sp + ".meta.tsv")) as mf:
            for m in csv.DictReader(mf, delimiter="\t"):
                per[m["sample"]].append(float(m["hcov"]) * float(m["mean_vcov_nonzero"]))
        depths[k].extend(sum(v) / len(v) for v in per.values())
print("samples in the MSA: species, genes removed by the coverage gate / genes in, median sample depth")
for k in sorted(by_n):
    s, rem, tot = by_n[k]
    dd = sorted(depths[k])
    print(f"  {k if k < 5 else '5+'}: {s} species, {rem}/{tot} ({rem / tot:.0%}), median depth {dd[len(dd) // 2]:.2f}x")
