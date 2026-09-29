#!/usr/bin/env python3
"""settled.py - species of the records a read's other genes settled (ZR:i:1), right or wrong, and
for the wrong ones how the read's confident genes voted."""

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
    by_read = collections.defaultdict(list)
    for line in open(args.sam):
        if line.startswith("@"):
            continue
        r = line.rstrip("\n").split("\t")
        if int(r[1]) & 0x100:
            continue
        bases = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHPX=])", r[5]) if op in "MXD=")
        by_read[r[0]].append((int(r[4]), names.get(r[2].split("_")[0], "?"), bases, "ZR:i:1" in r[11:]))
    right = wrong = 0
    wrong_votes = collections.Counter()
    for read, segs in by_read.items():
        for mapq, species, bases, settled in segs:
            if not settled:
                continue
            if species == truth[read]:
                right += bases
            else:
                wrong += bases
                confident = collections.Counter(s for m, s, _, z in segs if m >= 4 and not z)
                wrong_votes[f"confident genes: {dict(confident)}, true {truth[read]}"] += bases
    print(f"settled bp: right {right}, wrong {wrong} ({wrong / max(1, right + wrong):.4f})")
    for key, n in wrong_votes.most_common(8):
        print(n, key)


if __name__ == "__main__":
    main()
