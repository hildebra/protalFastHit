#!/usr/bin/env python3
"""How many of a protal SAM's fragments are duplicates (PCR or optical), per taxon: a fragment is its best record
pair's (gene, position, strand, mate's gene, mate's position), or for an unpaired record its (gene, position,
strand); fragments sharing that key with an earlier one are duplicates. Reads only primary records with a reference.
Prints per taxon (with at least --min-fragments fragments) the fragments, the distinct keys, the duplicate share, the
genes hit and the largest number of copies of one key, and the sample's totals.

usage: zstdcat SAMPLE.sam.zst | duplicates.py [--min-fragments 100] [--taxa 758,709]
"""
import argparse
import collections
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-fragments", type=int, default=100)
    ap.add_argument("--taxa", default="")
    args = ap.parse_args()
    wanted = set(args.taxa.split(",")) if args.taxa else None

    keys = collections.defaultdict(collections.Counter)  # taxon -> key -> copies
    genes = collections.defaultdict(set)
    seen = set()
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        f = line.split("\t", 9)
        flag = int(f[1])
        if flag & 0x904 or f[2] == "*":  # unmapped, secondary, supplementary
            continue
        if flag & 0x1 and flag & 0x80 and not (flag & 0x8):
            continue  # a pair's second record: counted with its first
        taxon, gene = f[2].split("_", 1)
        if wanted is not None and taxon not in wanted:
            continue
        name = f[0]
        if name in seen:
            continue
        seen.add(name)
        strand = "-" if flag & 0x10 else "+"
        mate = (f[6], f[7]) if flag & 0x1 and not (flag & 0x8) else ("", "")
        keys[taxon][(gene, f[3], strand) + mate] += 1
        genes[taxon].add(gene)

    total = distinct = 0
    print("taxon\tfragments\tdistinct\tduplicate_share\tgenes\tmost_copies")
    for taxon, counter in sorted(keys.items(), key=lambda kv: -sum(kv[1].values())):
        n = sum(counter.values())
        total += n
        distinct += len(counter)
        if n < args.min_fragments:
            continue
        print(f"{taxon}\t{n}\t{len(counter)}\t{1 - len(counter) / n:.3f}\t{len(genes[taxon])}\t{max(counter.values())}")
    if total:
        print(f"all\t{total}\t{distinct}\t{1 - distinct / total:.3f}")


if __name__ == "__main__":
    sys.exit(main())
