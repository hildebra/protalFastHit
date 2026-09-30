#!/usr/bin/env python3
"""crossmap.py - how many reads on each species' genes come from another species, and at what identity.

Usage: crossmap.py <protal_out_dir> <sample> [<sample> ...]
Read names carry the source genome (ART). Identity = M / (M + X + I + D) from the CIGAR (protal writes M
only for exact matches). Also reports whether the source species' own representative carries the gene
(a strain gene missing from its representative has nowhere else to map).
"""
import gzip, os, re, sys, pickle
from collections import defaultdict
import numpy as np

ACC = os.path.expanduser("~/audit5/accuracy")
CIG = re.compile(r"(\d+)([MIDNSHPX=])")
T = pickle.load(open(os.path.join(ACC, "truth.pkl"), "rb"))
GEN = T["meta"]["genomes"]
REPGENES = defaultdict(set)
for (sp, gid) in T["genes"]:
    REPGENES[sp].add(gid)
TAX = {}
for line in open(os.path.expanduser("~/audit5/world/protal_db/genome2tiid.tsv")):
    f = line.rstrip("\n").split("\t")
    TAX[f[1]] = f[3].split(";s__")[-1]


def main():
    out = sys.argv[1]
    ident = defaultdict(list)
    counts = defaultdict(int)
    for sample in sys.argv[2:]:
        with gzip.open(os.path.join(out, "alignments", sample + ".sam.gz"), "rt") as fh:
            for line in fh:
                if line[0] == "@":
                    continue
                f = line.split("\t", 7)
                if int(f[1]) & 4:
                    continue
                tax, gid = f[2].split("_")
                tsp = TAX[tax]
                src = f[0].split("_contig")[0]
                ssp = GEN[src]["species"]
                ops = defaultdict(int)
                for n, op in CIG.findall(f[5]):
                    ops[op] += int(n)
                idt = ops["M"] / max(1, ops["M"] + ops["X"] + ops["I"] + ops["D"])
                if ssp == tsp:
                    k = "own_species"
                else:
                    k = "other_species(src_rep_has_gene)" if int(gid) in REPGENES[ssp] else "other_species(src_rep_lacks_gene)"
                counts[k] += 1
                ident[k].append(idt)
    tot = sum(counts.values())
    for k in sorted(counts):
        a = np.array(ident[k])
        print(f"{k}: reads {counts[k]} ({counts[k] / tot:.4f}); identity median {np.median(a):.3f}, "
              f"5%/95% {np.percentile(a, 5):.3f}/{np.percentile(a, 95):.3f}; share <0.96: {(a < 0.96).mean():.3f}")


if __name__ == "__main__":
    main()
