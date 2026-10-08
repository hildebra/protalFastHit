#!/usr/bin/env python3
"""Depth over the simulated depth by gene length, for the species called and present in db_all (run in the validation
folder of run_validation.sh): results/depth_by_gene_length.txt."""
import csv, glob, collections, statistics
truth = collections.defaultdict(lambda: collections.defaultdict(float))
with open("sims/manifest.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        truth[row["sample"]][row["taxonomy"].split(";")[-1]] += float(row["vertical_coverage"])
bins = collections.defaultdict(list)
missing = collections.Counter()
for sample in truth:
    log = glob.glob(f"run_db_all/**/{sample}.profile.genes.log", recursive=True)[0]
    comp = glob.glob(f"run_db_all/**/{sample}.profile.composition", recursive=True)[0]
    share = float(next(csv.DictReader(open(comp), delimiter="\t"))["FragmentBaseShare"])
    with open(log) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row["Lineage"].split(";")[-1]
            if row["Predicted"] != "1" or name not in truth[sample]:
                continue
            length = int(row["GeneRefLength"])
            ratio = float(row["VCov"]) / (truth[sample][name] * share)
            b = min(length // 300, 6) * 300
            bins[b].append(ratio)
for b in sorted(bins):
    print(f"gene length {b}-{b+299}: median depth ratio {statistics.median(bins[b]):.4f} n={len(bins[b])}")
