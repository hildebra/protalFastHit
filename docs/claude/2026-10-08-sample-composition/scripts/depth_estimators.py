#!/usr/bin/env python3
"""Protal depth, the median and the length-weighted mean of the genes depths against the truth, by the median reads per gene
(db_all run of run_validation.sh; run in the validation folder): results/depth_estimators.txt."""
import csv, glob, collections, statistics
truth = collections.defaultdict(lambda: collections.defaultdict(float))
with open("sims/manifest.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        truth[row["sample"]][row["taxonomy"].split(";")[-1]] += float(row["vertical_coverage"])
rows = []
for sample in truth:
    comp = next(csv.DictReader(open(glob.glob(f"run_db_all/**/{sample}.profile.composition", recursive=True)[0]), delimiter="\t"))
    share = float(comp["FragmentBaseShare"])
    genes = collections.defaultdict(list)
    for row in csv.DictReader(open(glob.glob(f"run_db_all/**/{sample}.profile.genes.log", recursive=True)[0]), delimiter="\t"):
        name = row["Lineage"].split(";")[-1]
        if row["Predicted"] == "1" and name in truth[sample]:
            genes[name].append((float(row["VCov"]), int(row["GeneRefLength"]), int(row["TotalReads"])))
    log = {r["Name"]: float(r["VCov"]) for r in csv.DictReader(open(glob.glob(f"run_db_all/**/{sample}.profile.log", recursive=True)[0]), delimiter="\t")}
    for name, g in genes.items():
        true = truth[sample][name] * share
        med = statistics.median(v for v, _, _ in g)
        mean = sum(v * l for v, l, _ in g) / sum(l for _, l, _ in g)
        reads = statistics.median(n for _, _, n in g)
        rows.append((reads, med / true, mean / true, log[name] / true, len(g)))
def spread(v):
    v = sorted(v); p = lambda q: v[min(len(v) - 1, int(q * len(v)))]
    return f"{statistics.median(v):.4f} ({p(0.1):.3f}-{p(0.9):.3f}, n={len(v)})"
for label, lo, hi in (("median reads per gene < 10", 0, 10), ("10-30", 10, 30), (">= 30", 30, 10**9), ("all", 0, 10**9)):
    sel = [r for r in rows if lo <= r[0] < hi]
    if sel:
        print(f"{label}: protal depth {spread([r[3] for r in sel])}; median of genes {spread([r[1] for r in sel])}; length-weighted mean {spread([r[2] for r in sel])}")
