#!/usr/bin/env python3
"""Per database: how many genes the strain MSA's gene choice (SelectGenesForTaxon) leaves out.

usage: gene_uniques.py DIR...   (unpacked databases: reference.map, unique_kmers.tsv)
A gene is hittable with short + long unique k-mers > 0. The long-unique filter applies to a species
when 90% or more of its genes have long unique k-mers, and then drops its genes without any.
"""
import os
import sys
from collections import defaultdict


def main():
    for d in sys.argv[1:]:
        genes = defaultdict(set)
        with open(os.path.join(d, "reference.map")) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) >= 2 and f[0].isdigit():
                    genes[int(f[0])].add(int(f[1]))
        su, lu = {}, {}
        with open(os.path.join(d, "unique_kmers.tsv")) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                key = (int(f[0]), int(f[1]))
                su[key], lu[key] = int(f[2]), int(f[4])
        n_species = len(genes)
        n_genes = sum(len(g) for g in genes.values())
        no_unique = with_long = excluded = applies = 0
        excluded_frac = []
        for taxid, gs in genes.items():
            nu = sum(1 for g in gs if su.get((taxid, g), 0) + lu.get((taxid, g), 0) == 0)
            wl = sum(1 for g in gs if lu.get((taxid, g), 0) > 0)
            no_unique += nu
            with_long += wl
            ex = 0
            if wl * 10 >= len(gs) * 9:
                applies += 1
                ex = len(gs) - wl
            excluded += ex
            excluded_frac.append(ex / len(gs))
        excluded_frac.sort()
        q = lambda p: excluded_frac[min(len(excluded_frac) - 1, int(p * (len(excluded_frac) - 1)))]
        print(f"== {d}: {n_species} species, {n_genes} genes")
        print(f"  genes without any unique k-mer: {no_unique} ({no_unique / n_genes:.1%}); with long uniques: {with_long / n_genes:.1%}")
        print(f"  long-unique filter applies to {applies} species; genes it drops: {excluded} ({excluded / n_genes:.1%})")
        print(f"  dropped share per species: median {q(0.5):.1%}, p90 {q(0.9):.1%}, max {q(1.0):.1%}; "
              f"species losing >10%: {sum(f > 0.1 for f in excluded_frac)}, >25%: {sum(f > 0.25 for f in excluded_frac)}")


if __name__ == "__main__":
    main()
