#!/usr/bin/env python3
"""estimate_gene_rates.py DB_DIR MARKERS_TSV OUT_PREFIX [TRUE_RATES_TSV]

How fast each marker diverges within a species, relative to the genome, from what a database holds: every
genome's copy of each gene (full_reference.fna, >taxid_geneid) against the species representative's copy
(reference.fna). For each species and gene, the median k-mer (Mash, k = 12) distance of the other genomes'
copies to the representative's; divided by that species' median over its genes, so that how far a species'
strains are from its representative cancels out; the factor of a gene = the median of that ratio over the
species, scaled to a mean of 1. A gene's category (the simulator's marker_category of its name) gets the median
factor of its genes. Writes OUT_PREFIX.gene.tsv and OUT_PREFIX.category.tsv (geneid, factor), and with the
simulator's simulation/gene_rates.tsv also OUT_PREFIX.true.tsv, and prints how the estimates correlate with it.
"""
import collections
import csv
import math
import os
import statistics
import sys

db, markers_tsv, prefix = sys.argv[1:4]
truth = sys.argv[4] if len(sys.argv) > 4 else None
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "..", "scripts", "mini_db"))
from simulate_gtdb_release import marker_category  # noqa: E402

K = 12


def records(path):
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
            elif header:
                yield header, line.strip().upper()
                header = None


def kmers(seq):
    return {seq[i:i + K] for i in range(len(seq) - K + 1)}


def mash(a, b):
    union = len(a | b)
    j = len(a & b) / union if union else 0
    return 1.0 if j <= 0 else -math.log(2 * j / (1 + j)) / K


rep = {h: s for h, s in records(os.path.join(db, "reference.fna"))}
copies = collections.defaultdict(list)
for h, s in records(os.path.join(db, "full_reference.fna")):
    copies[h].append(s)
per_species = collections.defaultdict(dict)  # taxid -> geneid -> median distance of the other copies
for h, seqs in copies.items():
    if h not in rep:
        continue
    taxid, gid = h.split("_")
    r = kmers(rep[h])
    d = [mash(r, kmers(s)) for s in seqs if s != rep[h]]
    if d:
        per_species[taxid][int(gid)] = statistics.median(d)
ratios = collections.defaultdict(list)
for taxid, genes in per_species.items():
    if len(genes) < 10:
        continue
    typical = statistics.median(genes.values())
    if typical <= 0:
        continue
    for gid, d in genes.items():
        ratios[gid].append(d / typical)
factor = {gid: statistics.median(r) for gid, r in ratios.items() if len(r) >= 5}
mean = statistics.mean(factor.values())
factor = {g: f / mean for g, f in factor.items()}

name_of_marker = {}
with open(markers_tsv) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and not line.startswith(("#", "set\t")):
            name_of_marker[f[1]] = f[2]
marker_of = {}
with open(os.path.join(db, "gene2geneid.tsv")) as fh:
    for line in fh:
        marker, gid = line.rstrip("\n").split("\t")
        marker_of[int(gid)] = marker
category = {g: marker_category(name_of_marker.get(marker_of.get(g, ""), "")) for g in factor}
by_category = collections.defaultdict(list)
for g, f in factor.items():
    by_category[category[g]].append(f)
category_factor = {c: statistics.median(v) for c, v in by_category.items()}
with open(prefix + ".gene.tsv", "w") as fh:
    fh.writelines(f"{g}\t{f:.4f}\n" for g, f in sorted(factor.items()))
with open(prefix + ".category.tsv", "w") as fh:
    fh.writelines(f"{g}\t{category_factor[category[g]]:.4f}\n" for g in sorted(factor))
print(f"{len(factor)} genes from {len(per_species)} species; by category: " +
      ", ".join(f"{c} {category_factor[c]:.2f} ({len(v)} genes)" for c, v in sorted(by_category.items())))

if truth:
    true_rate = {}
    with open(truth) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            true_rate[row["marker"]] = float(row["rate"])
    gid_true = {g: true_rate[marker_of[g]] for g in factor if marker_of.get(g) in true_rate}
    with open(prefix + ".true.tsv", "w") as fh:
        fh.writelines(f"{g}\t{r:.4f}\n" for g, r in sorted(gid_true.items()))
    genes = sorted(gid_true)

    def corr(xs, ys):
        mx, my = statistics.mean(xs), statistics.mean(ys)
        num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        return num / math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))

    print(f"correlation with the true rates: per gene {corr([factor[g] for g in genes], [gid_true[g] for g in genes]):.3f}, "
          f"per category {corr([category_factor[category[g]] for g in genes], [gid_true[g] for g in genes]):.3f}")
