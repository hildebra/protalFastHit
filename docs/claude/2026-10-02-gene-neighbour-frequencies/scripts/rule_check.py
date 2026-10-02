#!/usr/bin/env python3
"""How rules for judging two genes next to each other on a read fare against the genomes' own gene order: for the
genes of each genome in a database's gene_positions.tsv, the verdicts the clade table (gene_neighbours.tsv) gives
- true pairs: two genes next to each other in the genome (within 3 kb), as a read across them shows them;
- skip pairs: a gene and the gene after its neighbour (within 3 kb), as a read whose middle gene was not aligned;
- far pairs: a gene and one of its genome's genes 10 kb or more away, as a chimeric read or a misplaced mate would.
A species' clades are its lineage's (internal_taxonomy.dmp). Rules (thresholds 0.2 / 0.05 on the share):
- nearest: the nearest clade with data; there, with the end informative in fewer than 3 species, a pairing seen is
  expected and one not seen is left to the next clade up (protal before the sparse fallback);
- fallback: the nearest clade with the end informative in 5 species or more decides; if none, the nearest with data
  says what it saw (protal's kMinInformative);
- smoothedM: the share of each clade pulled towards its parent's by M pseudo-species, from the top down
  ((species + M * parent share) / (informative + M)), at the nearest clade with data (M 2, 3, 5, 10).
Split by the size of the species' family (informative species at that gene end).

usage: rule_check.py <folder with gene_neighbours.tsv, gene_positions.tsv, internal_taxonomy.dmp> [--genomes N]
           [--positions FILE --held_out]
--held_out takes the genomes of --positions (a database with more species: the training database's full one) whose
species the folder's gene_positions.tsv lacks: species whose gene order the table does not know. --relatives judges
such a species' true pairs on the lineage of a species of the folder in its genus (else family), as the reads of a
relative the database lacks are, beside that species' own true pairs on its own lineage.
"""

import argparse
import collections
import os
import random

EXPECTED, UNLIKELY = 0.2, 0.05
MAX_GAP = 3000


def read_table(path):
    rules = collections.defaultdict(dict)  # (clade, gene, end) -> {(partner, partner end): species}
    informative = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if line.startswith("#") or f[0] == "clade":
                continue
            clade, gene, end, partner, partner_end, species, inf = map(int, f[:7])
            rules[(clade, gene, end)][(partner, partner_end)] = species
            informative[(clade, gene, end)] = inf
    return rules, informative


def verdict(share):
    return "expected" if share >= EXPECTED else "unlikely" if share <= UNLIKELY else "rare"


def judge(rule, chain, rules, informative, gene, end, partner, partner_end):
    """The verdict on `gene`'s `end` facing `partner`'s `partner_end` for a species of clades `chain` (nearest first)."""
    key = (partner, partner_end)
    if rule == "nearest":
        for clade in chain:
            k = (clade, gene, end)
            if k not in informative:
                continue
            seen = rules[k].get(key, 0)
            if informative[k] < 3:
                if seen:
                    return "expected"
                continue
            return verdict(seen / informative[k])
        return "unknown"
    if rule == "fallback":
        sparse = None
        for clade in chain:
            k = (clade, gene, end)
            if k not in informative:
                continue
            if informative[k] >= 5:
                return verdict(rules[k].get(key, 0) / informative[k])
            sparse = sparse or k
        if sparse is None:
            return "unknown"
        return "expected" if rules[sparse].get(key, 0) else "unknown"
    if rule.startswith("smoothed"):
        m = int(rule[len("smoothed"):])
        present = [c for c in chain if (c, gene, end) in informative]
        if not present:
            return "unknown"
        share = None
        for clade in reversed(present):  # from the top
            k = (clade, gene, end)
            seen = rules[k].get(key, 0)
            share = seen / informative[k] if share is None else (seen + m * share) / (informative[k] + m)
        return verdict(share)
    raise ValueError(rule)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("db")
    ap.add_argument("--genomes", type=int, default=600, help="genomes to sample (default 600)")
    ap.add_argument("--positions", help="the gene positions of the genomes to judge (default: the folder's)")
    ap.add_argument("--held_out", action="store_true", help="only the species that the folder's positions lack")
    ap.add_argument("--relatives", action="store_true", help="with --held_out: on a relative's lineage (above)")
    args = ap.parse_args()
    rules, informative = read_table(os.path.join(args.db, "gene_neighbours.tsv"))
    clades = {c for c, _, _ in informative}
    nodes = {}
    with open(os.path.join(args.db, "internal_taxonomy.dmp")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            nodes[int(f[0])] = (int(f[1]), f[4])

    def chain(taxid):
        out, t = [], taxid
        while True:
            if t in clades:
                out.append(t)
            if nodes[t][0] == t:
                return out
            t = nodes[t][0]

    def family(taxid):
        t = taxid
        while nodes[t][1] != "family" and nodes[t][0] != t:
            t = nodes[t][0]
        return t

    def read_genomes(path, skip=frozenset()):
        out = collections.defaultdict(list)
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if line.startswith("#") or f[0] == "accession" or f[1] in skip:
                    continue
                out[(f[0], int(f[1]))].append((f[2], int(f[6]), int(f[7]), int(f[5]), f[8]))
        return out

    own = read_genomes(os.path.join(args.db, "gene_positions.tsv"))
    known = {str(t) for _, t in own}
    genomes = read_genomes(args.positions, known) if args.held_out else own
    rng = random.Random(1)
    sample = rng.sample(sorted(genomes), min(args.genomes, len(genomes)))
    jobs = [(genomes[g], g[1], "") for g in sample]
    if args.relatives:
        by_parent = collections.defaultdict(list)
        for acc, t in sorted(own):
            by_parent[nodes[t][0]].append((acc, t))
        jobs = []
        for acc, taxid in sample:
            genus = nodes[taxid][0]
            near = by_parent.get(genus) or [g for p, gs in by_parent.items() if family(p) == family(genus) for g in gs]
            if not near:
                continue
            relative = rng.choice(near)
            jobs.append((genomes[(acc, taxid)], relative[1], "relative's "))
            jobs.append((own[relative], relative[1], "own "))
    rules_judged = ("nearest", "fallback", "smoothed2", "smoothed3", "smoothed5", "smoothed10")
    counts = collections.defaultdict(collections.Counter)  # (rule, kind, family size class) -> verdicts
    for positions, taxid, label in jobs:
        lineage = chain(taxid)
        fam = family(taxid)
        by_contig = collections.defaultdict(list)
        for contig, start, end, gene, strand in positions:
            by_contig[contig].append((start, end, gene, strand))
        everything = sorted(g for genes in by_contig.values() for g in genes)
        for genes in by_contig.values():
            genes.sort()
            for i, (start, end, gene, strand) in enumerate(genes):
                right = 3 if strand == "+" else 5
                pairs = []
                for kind, j in (("true", i + 1), ("skip", i + 2)):
                    if j < len(genes) and genes[j][0] - end <= MAX_GAP:
                        pairs.append((kind, genes[j]))
                far = [g for g in everything if g[2] != gene and (g[0] - end > 10000 or start - g[1] > 10000)]
                if far:
                    pairs.append(("far", rng.choice(far)))
                for kind, (_, _, partner, partner_strand) in pairs:
                    kind = label + kind
                    partner_end = 5 if partner_strand == "+" else 3
                    n = informative.get((fam, gene, right), 0)
                    size = "family >= 5" if n >= 5 else "family < 5"
                    for rule in rules_judged:
                        v = judge(rule, lineage, rules, informative, gene, right, partner, partner_end)
                        counts[(rule, kind, size)][v] += 1
                        counts[(rule, kind, "all")][v] += 1
    print(f"{len(sample)} genomes of {len({t for _, t in sample})} species" +
          (f"; {len(jobs) // 2} with a relative in the folder" if args.relatives else ""))
    print("| pairs | family | rule | pairs judged | expected | rare | unlikely | unknown |")
    print("|---|---|---|---|---|---|---|---|")
    kinds = [k for k in ("true", "skip", "far", "own true", "relative's true", "own skip", "relative's skip")
             if any(key[1] == k for key in counts)]
    for kind in kinds:
        for size in ("all", "family >= 5", "family < 5"):
            for rule in rules_judged:
                c = counts[(rule, kind, size)]
                total = sum(c.values())
                if total:
                    print(f"| {kind} | {size} | {rule} | {total} | " +
                          " | ".join(f"{100 * c[v] / total:.1f}%" for v in ("expected", "rare", "unlikely", "unknown")) + " |")


if __name__ == "__main__":
    main()
