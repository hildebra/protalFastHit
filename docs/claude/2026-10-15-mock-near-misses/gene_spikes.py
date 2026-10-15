#!/usr/bin/env python3
"""Taxa whose reads pile on one or two genes far beyond the rest of their genes: per taxon of a protal SAM, the
primary records per gene, the median over the genes hit, the top gene's records and their share, and how many genes
exceed a Poisson outlier bound around the median (P < 1e-6 per gene). Prints the taxa with at least --min-records
records whose top gene holds more than --min-share of them.

usage: zstdcat SAMPLE.sam.zst | gene_spikes.py [--min-records 50] [--min-share 0.3] [--names genome2tiid.tsv]
"""
import argparse
import collections
import math
import statistics
import sys


def poisson_upper(mean, p=1e-6):
    """The smallest k with P(X >= k) < p for X ~ Poisson(mean)."""
    k, term = 0, math.exp(-mean)
    cdf = term
    while 1 - cdf >= p and k < 10_000:
        k += 1
        term *= mean / k
        cdf += term
    return k + 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-records", type=int, default=50)
    ap.add_argument("--min-share", type=float, default=0.3)
    ap.add_argument("--names")
    args = ap.parse_args()
    names = {}
    if args.names:
        with open(args.names) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) >= 4:
                    names[f[1]] = f[3].split(";")[-1]

    per_gene = collections.defaultdict(collections.Counter)
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.split("\t", 4)
        if int(f[1]) & 0x904 or f[2] == "*":
            continue
        taxon, gene = f[2].split("_", 1)
        per_gene[taxon][gene] += 1

    print("taxon\tname\trecords\tgenes\tmedian\ttop_gene\ttop_records\ttop_share\toutlier_genes\trecords_without_outliers")
    for taxon, genes in sorted(per_gene.items(), key=lambda kv: -sum(kv[1].values())):
        n = sum(genes.values())
        if n < args.min_records:
            continue
        top_gene, top = genes.most_common(1)[0]
        if top / n < args.min_share:
            continue
        median = statistics.median(genes.values())
        bound = poisson_upper(max(median, 1.0))
        outliers = [g for g, c in genes.items() if c >= bound]
        rest = n - sum(genes[g] for g in outliers)
        print(f"{taxon}\t{names.get(taxon, '')}\t{n}\t{len(genes)}\t{median:g}\t{top_gene}\t{top}\t{top / n:.3f}\t"
              f"{len(outliers)}\t{rest}")


if __name__ == "__main__":
    sys.exit(main())
