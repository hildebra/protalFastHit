#!/usr/bin/env python3
"""foreign_world_others.py OUT - the foreign genes of foreign_world.sh's kept runs other than the donor's gene A: on a
species of the sample (the donor or a background species: a false flag) or on a taxon that is not in it (one holding
reads of the recipient, an unknown species, or of a relative), and whether the taxon is called.
"""
import collections
import csv
import os
import sys

out = sys.argv[1]
reps = {}
with open(os.path.join(out, "full_species", "internal_taxonomy.dmp")) as fh:
    next(fh)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) > 6 and f[4] == "species" and f[6]:
            reps[f[6]] = f[0]
triples = {r["sample"].removeprefix("hgt"): r for r in csv.DictReader(open(os.path.join(out, "world", "triples.tsv")), delimiter="\t")}
present = collections.defaultdict(set)  # sample -> DB taxids of its genomes
for r in csv.DictReader(open(os.path.join(out, "world", "samples.tsv")), delimiter="\t"):
    name = os.path.basename(r["genome"])[:-len(".fna")]
    i = "".join(c for c in r["sample"] if c.isdigit())
    if name.startswith("donor"):
        present[r["sample"]].add(triples[i]["t_taxid"])
    elif name.startswith("GC"):
        present[r["sample"]].add(reps.get(name, "?"))
by_kind = collections.Counter()
print("| reads | foreign genes on species of the sample (called / not) | on taxa not in it (called / not) | judged genes of present species |")
print("|---|---|---|---|")
for kind in ("pe", "pb", "ont"):
    counts = collections.Counter()
    judged_present = 0
    for sample in sorted(present):
        path = os.path.join(out, "runs", f"{kind}.keep", f"{sample}.profile.genes.log")
        if not os.path.exists(path):
            continue
        i = "".join(c for c in sample if c.isdigit())
        for r in csv.DictReader(open(path), delimiter="\t"):
            mine = r["TaxID"] in present[sample]
            judged_present += mine and int(r["LinksJudged"]) >= 4
            if r["Foreign"] != "1":
                continue
            if sample.startswith("hgt") and r["TaxID"] == triples[i]["t_taxid"] and r["GeneID"] == triples[i]["gene"]:
                continue
            counts[(mine, r["Predicted"] == "1")] += 1
            if mine:
                by_kind[(kind, sample[:3])] += 1
    print(f"| {kind} | {counts[(True, True)]} / {counts[(True, False)]} | {counts[(False, True)]} / {counts[(False, False)]} | "
          f"{judged_present} |")
print()
print("foreign genes on species of the sample, by sample kind:", dict(sorted(by_kind.items())))
