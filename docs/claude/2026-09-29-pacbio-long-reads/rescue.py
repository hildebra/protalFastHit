#!/usr/bin/env python3
"""rescue.py - how much of a native run's ambiguous genes (MAPQ < 4) a read's other genes could settle.

For each read, its segments (primary and supplementary records) with MAPQ >= 4 name species; if they
all name one species, the read's ambiguous segments could be given to it. Reports the ambiguous
bases on such reads and how often that species is the read's true one.
"""

import argparse
import collections
import re


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sam", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--taxonomy", required=True)
    args = ap.parse_args()
    names = {}
    for line in open(args.taxonomy):
        f = line.rstrip("\n").split("\t")
        if f[0].isdigit():
            names[f[0]] = f[3]
    truth = {}
    for line in open(args.truth):
        f = line.rstrip("\n").split("\t")
        if f[0] != "read":
            truth[f[0]] = f[2]
    segments = collections.defaultdict(list)
    for line in open(args.sam):
        if line.startswith("@"):
            continue
        r = line.split("\t", 10)
        if int(r[1]) & 0x100:
            continue
        bases = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHPX=])", r[5]) if op in "MXD=")
        segments[r[0]].append((int(r[4]), names.get(r[2].split("_")[0], "?"), bases))
    ambiguous = settled = settled_right = no_confident = mixed = 0
    for read, segs in segments.items():
        confident = {species for mapq, species, _ in segs if mapq >= 4}
        amb = sum(b for mapq, _, b in segs if mapq < 4)
        ambiguous += amb
        if not confident:
            no_confident += amb
        elif len(confident) > 1:
            mixed += amb
        else:
            settled += amb
            settled_right += amb if confident == {truth[read]} else 0
    print(f"ambiguous bp {ambiguous}; on reads whose confident genes name one species: {settled} "
          f"({settled / max(1, ambiguous):.3f}), that species right for {settled_right / max(1, settled):.4f}; "
          f"reads without confident genes: {no_confident}; with several species: {mixed}")


if __name__ == "__main__":
    main()
