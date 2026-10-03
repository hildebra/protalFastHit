#!/usr/bin/env python3
"""The suspect gene copies of the r226 build (gene_incongruence.tsv of v7): by the rank the pair shares, by gene,
by species; how many species lose how many of their marker genes."""
import csv, sys
from collections import Counter, defaultdict

path, taxonomy = sys.argv[1], sys.argv[2]
names = {}
with open(taxonomy, encoding="utf-8") as fh:  # internal_taxonomy.dmp: id, parent_id, external_id, name, rank, level, rep_genome (tabs)
    next(fh)
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 4:
            names[parts[0]] = parts[3]

suspects = set()
suspect_rank = Counter()
suspect_dist = Counter()
pair_rank = Counter()
near_identical = Counter()
genes_of_species = defaultdict(set)
partners = defaultdict(Counter)
identical_partner = defaultdict(int)
with open(path, encoding="utf-8") as fh:
    rd = csv.DictReader(fh, delimiter="\t")
    for r in rd:
        rank = r["shared_rank"]
        d = abs(float(r["distance"]))
        pair_rank[rank] += 1
        if d <= 0.001:
            near_identical[rank] += 1
        for side, other in (("a", "b"), ("b", "a")):
            if r[f"suspect_{side}"] == "1":
                key = (r["geneid"], r[f"taxid_{side}"])
                if key not in suspects:
                    suspects.add(key)
                    suspect_rank[rank] += 1
                    suspect_dist["0" if d <= 0.001 else "<=0.01" if d <= 0.01 else "<=0.02"] += 1
                    genes_of_species[r[f"taxid_{side}"]].add(r["geneid"])
                    if d <= 0.001:
                        identical_partner[r[f"taxid_{side}"]] += 1
                partners[r[f"taxid_{side}"]][r[f"taxid_{other}"]] += 1

print(f"near pairs across genera (within 0.05): {sum(pair_rank.values())}; by the rank shared: {dict(pair_rank.most_common())}")
print(f"of them identical (distance 0): {sum(near_identical.values())}: {dict(near_identical.most_common())}")
print(f"\nsuspect copies: {len(suspects)} of {len(genes_of_species)} species")
print(f"by the rank the flagged pair shares: {dict(suspect_rank.most_common())}")
print(f"by distance to the other genus's copy: {dict(suspect_dist.most_common())}")
per_species = Counter({s: len(g) for s, g in genes_of_species.items()})
order = ["1", "2", "3-5", "6-20", "21-60", ">60"]
hist = Counter()
for s, n in per_species.items():
    hist["1" if n == 1 else "2" if n == 2 else "3-5" if n <= 5 else "6-20" if n <= 20 else "21-60" if n <= 60 else ">60"] += 1
print(f"\nspecies by the number of their genes flagged: {[(k, hist[k]) for k in order if hist[k]]}")
print("\nthe 30 species with the most genes flagged (of 120 bac / 53 ar markers), their copies identical to the partner's, and the partner species:")
for s, n in per_species.most_common(30):
    top = partners[s].most_common(2)
    print(f"  {n:4d} ({identical_partner[s]:3d} identical)  {names.get(s, s):42s}  partners: " + "; ".join(f"{names.get(t, t)} ({c})" for t, c in top))
gene_counts = Counter(g for g, _ in suspects)
print(f"\ngenes by flagged copies: top 10 {gene_counts.most_common(10)}; genes with none flagged: {168 - len(gene_counts)}")
