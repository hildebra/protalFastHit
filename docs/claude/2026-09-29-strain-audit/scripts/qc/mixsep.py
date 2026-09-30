#!/usr/bin/env python3
"""Per sample row of a run: truly mixed (>=2 strains each >= 1x) or not, depth, the count qcmsa
fences (genes with MRate2 > 0) and the pooled multi-allelic site rate (sum multi_allelic / sum counts_vcov2)."""
import sys, os, csv, glob
from collections import defaultdict

run = sys.argv[1]
man = defaultdict(list)
for r in csv.DictReader(open(os.path.join(run, "sim", "manifest.tsv")), delimiter="\t"):
    man[("s__" + r["species"].replace(" ", "_"), r["sample"])].append(float(r["vertical_coverage"]))
rows = []
for meta in glob.glob(os.path.join(run, "protal", "strains", "*.meta.tsv")):
    sp = os.path.basename(meta)[:-9]
    acc = defaultdict(lambda: [0, 0, 0, 0.0, 0])
    for r in csv.DictReader(open(meta), delimiter="\t"):
        a = acc[r["sample"]]
        a[0] += float(r["multi_rate_vcov2"]) > 0
        a[1] += int(r["multi_allelic"])
        a[2] += int(r["counts_vcov2"])
        a[3] += float(r["vertical_coverage"]); a[4] += 1
    summ = os.path.join(run, "protal", "strains", sp + ".qcmsa_summary.tsv")
    removed = set()
    if os.path.exists(summ):
        removed = {l.split("\t")[1] for l in open(summ) if l.startswith("sample_filtered")}
    for s, a in acc.items():
        covs = man.get((sp, s), [])
        mixed = sum(c >= 1 for c in covs) >= 2
        rows.append((sp, s, mixed, a[3] / a[4], a[0], a[1] / a[2] if a[2] else 0, s in removed))
for mixed in (True, False):
    sel = [r for r in rows if r[2] == mixed]
    rates = sorted(r[5] for r in sel)
    if not rates: continue
    print(f"{'mixed' if mixed else 'single-strain'}: n={len(sel)} removed by qcmsa={sum(r[6] for r in sel)} "
          f"multi-site rate min {rates[0]:.4%} median {rates[len(rates)//2]:.4%} max {rates[-1]:.4%}")
for r in sorted(rows, key=lambda r: -r[5])[:16]:
    print(f"  {r[0]:22} {r[1]:10} mixed={r[2]!s:5} depth={r[3]:6.1f} multi_genes={r[4]:4} multi_rate={r[5]:.4%} removed={r[6]}")
