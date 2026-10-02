#!/usr/bin/env python3
"""Which rank decides on the gene ends of a database's species, as protal's gene_neighbours::Table::Deciding picks it:
the nearest clade with the end informative in MIN species or more, else the nearest with any data ("sparse"), for
MIN 1 (the nearest clade with data, as before the sparse fallback), 3 and 5 (kMinInformative).

usage: deciding_rank.py <converted folder>   (gene_neighbours.tsv, internal_taxonomy.dmp, reference.map)
"""

import collections
import os
import sys


def main():
    db = sys.argv[1]
    nodes = {}
    with open(os.path.join(db, "internal_taxonomy.dmp")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            nodes[int(f[0])] = (int(f[1]), f[4], f[3])
    informative = {}  # (clade, gene, end) -> informative species
    with open(os.path.join(db, "gene_neighbours.tsv")) as fh:
        for line in fh:
            f = line.split("\t")
            if line.startswith("#") or f[0] == "clade":
                continue
            informative[(int(f[0]), int(f[1]), int(f[2]))] = int(f[6])
    clades = {c for c, _, _ in informative}
    genes = collections.defaultdict(list)
    with open(os.path.join(db, "reference.map")) as fh:
        for line in fh:
            f = line.split()
            genes[int(f[0])].append(int(f[1]))

    def chain(taxid):
        out, t = [], taxid
        while True:
            parent = nodes[t][0]
            if t in clades:
                out.append(t)
            if parent == t:
                return out
            t = parent

    def domain(taxid):
        t = taxid
        while nodes[t][1] != "domain" and nodes[t][0] != t:
            t = nodes[t][0]
        return t

    for minimum in (1, 3, 5):
        by_domain = collections.defaultdict(collections.Counter)
        for taxid, ids in genes.items():
            lineage = chain(taxid)
            name = nodes[domain(taxid)][2]
            for gene in ids:
                for end in (5, 3):
                    decided, sparse = None, None
                    for clade in lineage:
                        n = informative.get((clade, gene, end))
                        if n is None:
                            continue
                        if n >= minimum:
                            decided = nodes[clade][1]
                            break
                        sparse = sparse or nodes[clade][1]
                    by_domain[name][decided or (f"sparse {sparse}" if sparse else "none")] += 1
        for name, counts in sorted(by_domain.items()):
            total = sum(counts.values())
            order = ["family", "order", "class", "phylum", "domain"]
            keys = order + sorted(k for k in counts if k not in order)
            print(f"min {minimum}, domain {name}: {total} gene ends: " +
                  ", ".join(f"{k} {100 * counts[k] / total:.1f}%" for k in keys if counts[k]))


if __name__ == "__main__":
    main()
