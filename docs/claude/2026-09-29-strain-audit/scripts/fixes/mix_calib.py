#!/usr/bin/env python3
"""Pooled multi-allelic rates of single-strain vs mixed MSA rows in the accuracy runs, for a label.

usage: mix_calib.py LABEL [FLOOR...]
The strains of a (sample, species) come from sim_<run>/manifest.tsv (two genomes: a mixture); the
rates (sum multi_allelic / sum counts_vcov2) from prot_<run>_LABEL/strains/<species>.meta.tsv.
For each FLOOR: the mixtures and single-strain rows above it.
"""
import csv
import os
import sys
from collections import defaultdict

A = os.path.expanduser("~/audit5/accuracy")
label = sys.argv[1]
floors = [float(f) for f in sys.argv[2:]] or [0.004, 0.002, 0.001, 0.0005]
single, mixed = [], []
for run in ("A", "B", "Cs", "Cl"):
    man = os.path.join(A, f"sim_{'C' if run.startswith('C') else run}", "manifest.tsv")
    strains = os.path.join(A, f"prot_{run}_{label}", "strains")
    if not os.path.exists(man) or not os.path.isdir(strains):
        continue
    depth = defaultdict(list)
    with open(man) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            depth[(r["sample"], r["species"])].append(float(r["vertical_coverage"]))
    for f in os.listdir(strains):
        if not f.endswith(".meta.tsv"):
            continue
        species = f[len("s__"):-len(".meta.tsv")].replace("_", " ")
        multi, pos, bad = defaultdict(int), defaultdict(int), defaultdict(int)
        with open(os.path.join(strains, f)) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                multi[r["sample"]] += int(r["multi_allelic"])
                pos[r["sample"]] += int(r["counts_vcov2"])
                bad[r["sample"]] += int(r["multi_allelic"]) > 0
        for s in multi:
            d = depth.get((s, species))
            if not d:
                continue
            rate = multi[s] / pos[s] if pos[s] else 0.0
            item = (rate, bad[s], run, s, species, sum(d), min(d) / sum(d))
            (mixed if len(d) > 1 else single).append(item)

single.sort(reverse=True)
mixed.sort(reverse=True)
print(f"== single-strain rows: {len(single)}; highest:")
for rate, bad, run, s, sp, d, _ in single[:8]:
    print(f"  {rate:.5f}  genes {bad:3d}  {run} {s} {sp} depth {d:.3g}")
print(f"== mixed rows: {len(mixed)}")
for rate, bad, run, s, sp, d, minor in mixed:
    print(f"  {rate:.5f}  genes {bad:3d}  {run} {s} {sp} depth {d:.3g} minor {minor:.2f}")
for floor in floors:
    m = [x for x in mixed if x[0] > floor]
    s_ = [x for x in single if x[0] > floor]
    print(f"floor {floor}: mixtures above {len(m)}/{len(mixed)} (minor >= 0.15: "
          f"{sum(x[6] >= 0.15 for x in m)}/{sum(x[6] >= 0.15 for x in mixed)}), single-strain above {len(s_)}/{len(single)}")
