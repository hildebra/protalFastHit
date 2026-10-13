#!/usr/bin/env python3
"""The long reads whose consensus taxon differs between two protal SAMs of the same sample: per read the taxon most of
its records name (those not tagged ZR:i:2, as ../2026-10-13-ecoli-hifi/lost_reads.py), and of the reads that changed,
how many moved from a species of the truth to one outside it, the other way, between two of the truth, or between two
outside it. Reads with records in one SAM only are counted apart.

usage: consensus_changes.py SAM_A SAM_B GENOME2TIID TRUTH_TAXIDS [READS]
  (TRUTH_TAXIDS: one taxid per line, ../2026-10-13-ecoli-hifi/truth_taxids.py; READS: list up to this many read names
  of each of the most frequent moves)
"""
import collections
import subprocess
import sys


def consensus(path):
    taxa = collections.defaultdict(collections.Counter)
    text = subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE, text=True).stdout
    for line in text:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if int(f[1]) & (4 | 256) or "ZR:i:2" in f[11:]:
            continue
        taxa[f[0]][f[2].split("_", 1)[0]] += 1
    return {read: c.most_common(1)[0][0] for read, c in taxa.items()}


def main():
    a, b = consensus(sys.argv[1]), consensus(sys.argv[2])
    names = {}
    with open(sys.argv[3]) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            names[f[1]] = f[3].rsplit("s__", 1)[-1]
    with open(sys.argv[4]) as fh:
        truth = {line.split()[0] for line in fh if line.strip()}
    show = int(sys.argv[5]) if len(sys.argv) > 5 else 0
    kinds = collections.Counter()
    moves = collections.Counter()
    examples = collections.defaultdict(list)
    for read in sorted(a.keys() & b.keys()):
        x, y = a[read], b[read]
        if x == y:
            kinds["same"] += 1
            continue
        kind = {(True, False): "truth -> other", (False, True): "other -> truth", (True, True): "truth -> truth",
                (False, False): "other -> other"}[(x in truth, y in truth)]
        kinds[kind] += 1
        move = f"{names.get(x, x)} -> {names.get(y, y)}"
        moves[move] += 1
        if len(examples[move]) < show:
            examples[move].append(read)
    print(f"# {sys.argv[1]} -> {sys.argv[2]}: reads in both {len(a.keys() & b.keys())}, only in the first "
          f"{len(a.keys() - b.keys())}, only in the second {len(b.keys() - a.keys())}")
    for kind in ("same", "other -> truth", "truth -> other", "truth -> truth", "other -> other"):
        print(f"{kind}\t{kinds[kind]}")
    print("# most frequent moves")
    for move, n in moves.most_common(10):
        print(f"{n}\t{move}" + (f"\t{' '.join(examples[move])}" if show else ""))


if __name__ == "__main__":
    main()
