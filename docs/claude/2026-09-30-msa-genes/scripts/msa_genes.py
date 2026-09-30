#!/usr/bin/env python3
"""Which of a species' genes reach its strain MSAs, and why the others do not.

usage: msa_genes.py RUN_DIR DB_DIR   (RUN_DIR: protal output with strains/ and a protal.log beside it;
DB_DIR: the unpacked database, reference.map and unique_kmers.tsv)
Per species with an MSA (strains/species.tsv): the genome's genes, those with reads in its MSA
samples (.meta.tsv), in the raw MSA (.raw.partition.txt) and after qcmsa (.partition.txt); the
missing genes by cause. Also: MSA genes by their unique k-mer count.
"""
import csv
import os
import re
import sys
from collections import Counter, defaultdict

run, db = sys.argv[1], sys.argv[2]
strains = os.path.join(run, "strains")
genes = defaultdict(set)
with open(os.path.join(db, "reference.map")) as fh:
    for line in fh:
        f = line.split("\t")
        genes[int(f[0])].add(int(f[1]))
uniq = {}
with open(os.path.join(db, "unique_kmers.tsv")) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        uniq[(int(f[0]), int(f[1]))] = (int(f[2]) + int(f[4]), int(f[8]))
log = open(os.path.join(os.path.dirname(run.rstrip("/")), "protal.log")).read() if os.path.exists(
    os.path.join(os.path.dirname(run.rstrip("/")), "protal.log")) else ""
long_unique_drop = {m.group(1): int(m.group(2)) for m in
                    re.finditer(r"^(\S+): (\d+) of \d+ genes have no long unique k-mers", log, re.M)}


def partition_genes(path):
    out = set()
    if os.path.exists(path):
        for line in open(path):
            m = re.match(r"DNA, gene(\d+) =", line)
            if m:
                out.add(int(m.group(1)))
    return out


totals = Counter()
per_species = []
by_uniques = defaultdict(lambda: [0, 0])  # unique k-mer bin -> [genes, in filtered MSA]
with open(os.path.join(strains, "species.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        sp, taxid = r["species"], int(r["taxid"])
        if r["raw_msa"] == "-":
            continue
        genome = genes[taxid]
        with_reads = set()
        with open(os.path.join(strains, sp + ".meta.tsv")) as mf:
            for m in csv.DictReader(mf, delimiter="\t"):
                with_reads.add(int(m["gene_id"]))
        raw = partition_genes(os.path.join(strains, sp + ".raw.partition.txt"))
        filt = partition_genes(os.path.join(strains, sp + ".partition.txt"))
        summary = {}
        spath = os.path.join(strains, sp + ".qcmsa_summary.tsv")
        if os.path.exists(spath):
            for line in open(spath):
                f = line.rstrip("\n").split("\t")
                if f[0] == "gene_filtered":
                    summary[int(f[1])] = "coverage" if "coverage" in f[3] else "multi-allelic"
        lu = long_unique_drop.get(sp, 0)
        row = Counter(genome=len(genome), reads=len(with_reads & genome), raw=len(raw), filt=len(filt),
                      no_reads=len(genome - with_reads), long_unique=lu,
                      no_depth=len(with_reads & genome) - lu - len(raw),
                      qc_coverage=sum(1 for g in raw - filt if summary.get(g) == "coverage"),
                      qc_multi=sum(1 for g in raw - filt if summary.get(g) == "multi-allelic"))
        row["qc_sites"] = len(raw - filt) - row["qc_coverage"] - row["qc_multi"]
        totals.update(row)
        per_species.append((len(filt) / len(genome), sp, int(r["samples"]), row))
        for g in genome:
            u, total = uniq.get((taxid, g), (0, 0))
            b = "0" if u == 0 else "1-9" if u < 10 else "10-99" if u < 100 else "100+"
            by_uniques[b][0] += 1
            by_uniques[b][1] += g in filt

n = len(per_species)
print(f"{n} species with an MSA; genes: {totals['genome']} in their genomes")
for key, label in (("reads", "with reads in the MSA's samples"), ("raw", "in the raw MSAs"), ("filt", "after qcmsa")):
    print(f"  {label}: {totals[key]} ({totals[key] / totals['genome']:.1%})")
print("  missing, by cause: " + ", ".join(f"{k} {totals[k]}" for k in
      ("no_reads", "long_unique", "no_depth", "qc_coverage", "qc_multi", "qc_sites")))
per_species.sort()
print("species with the fewest genes after qcmsa (share, samples, genome/reads/raw/filtered):")
for frac, sp, samples, row in per_species[:12]:
    print(f"  {frac:.0%} {sp} ({samples} samples): {row['genome']}/{row['reads']}/{row['raw']}/{row['filt']}"
          f"  no reads {row['no_reads']}, long-unique {row['long_unique']}, qc coverage {row['qc_coverage']}, qc multi {row['qc_multi']}")
fr = sorted(f for f, *_ in per_species)
print(f"share of genome genes after qcmsa per species: median {fr[n // 2]:.0%}, p10 {fr[n // 10]:.0%}, min {fr[0]:.0%}")
print("genes by unique k-mers (short + long): in the filtered MSAs / all")
for b in ("0", "1-9", "10-99", "100+"):
    g, f = by_uniques[b]
    if g:
        print(f"  {b:>6}: {f}/{g} ({f / g:.0%})")
