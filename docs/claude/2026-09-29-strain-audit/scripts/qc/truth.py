#!/usr/bin/env python3
"""Pairwise p-distances in an MSA vs the simulated truth (star tree: d(i,j) = d_i + d_j)."""
import sys, os, csv, itertools
sys.path.insert(0, os.path.dirname(__file__))
from consistency import read_fasta

WORLD = os.path.expanduser("~/audit5/world/gtdb_r226/simulation")

def load_truth(manifest):
    div = {r["accession"]: float(r["strain_divergence"]) for r in csv.DictReader(open(os.path.join(WORLD, "divergence.tsv")), delimiter="\t")}
    rep = {}
    for r in csv.DictReader(open(os.path.join(WORLD, "genomes.tsv")), delimiter="\t"):
        if r["gtdb_representative"] == "t":
            sp = "s__" + r["gtdb_taxonomy"].split(";s__")[1].replace(" ", "_")
            rep[sp] = r["accession"]
    samp = {}
    for r in csv.DictReader(open(manifest), delimiter="\t"):
        sp = "s__" + r["species"].replace(" ", "_")
        samp.setdefault(sp, {}).setdefault(r["sample"], []).append(r["genome"])
    return div, rep, samp

def genome_of(name, sp, rep, samp):
    if name.endswith("_reference"):
        return [rep[sp]]
    return samp.get(sp, {}).get(name, ["?"])

def true_dist(g1, g2, div):
    if len(g1) != 1 or len(g2) != 1 or "?" in g1 or "?" in g2:
        return float("nan")
    if g1 == g2:
        return 0.0
    return div[g1[0]] + div[g2[0]]

def pdist(a, b):
    n = d = 0
    for x, y in zip(a, b):
        if x in "ACGT" and y in "ACGT":
            n += 1
            d += x != y
    return (d / n if n else float("nan")), n

def main(msa, sp, manifest):
    div, rep, samp = load_truth(manifest)
    names, seqs = read_fasta(msa)
    print(f"{sp}: {msa}")
    for (i, a), (j, b) in itertools.combinations(enumerate(names), 2):
        g1, g2 = genome_of(a, sp, rep, samp), genome_of(b, sp, rep, samp)
        p, n = pdist(seqs[i], seqs[j])
        t = true_dist(g1, g2, div)
        print(f"  {a:>32} {b:>10}  p={p:.5f} (n={n:>6})  truth={t:.5f}  {'SAME' if g1 == g2 else ''}")

if __name__ == "__main__":
    main(*sys.argv[1:4])
