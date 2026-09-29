#!/usr/bin/env python3
"""diagnose.py - marker genes a native run misses that the 1 kb split run finds, and why.

For each read, the marker genes of its genome it overlaps (marker_positions.tsv, >= 200 bp) are
listed with how the native SAM covers them: found (a counted record overlapping it on the read),
or missed. Missed genes are summarised by where they lie on the read (cut by a read end or not)
and whether the split run found a record for the same read region.
"""

import argparse
import collections
import re


def read_table(path):
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip() and not l.startswith("#")]
    header = lines[0].split("\t")
    return [dict(zip(header, l.split("\t"))) for l in lines[1:]]


def native_intervals(sam, lengths):
    """Per read, the read intervals (forward coordinates) of its counted records."""
    out = collections.defaultdict(list)
    for line in open(sam):
        if line.startswith("@"):
            continue
        r = line.split("\t", 10)
        flag, mapq = int(r[1]), int(r[4])
        if flag & 0x100 or mapq < 4:
            continue
        ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHPX=])", r[5])]
        left = ops[0][0] if ops[0][1] == "H" else 0
        query = sum(n for n, op in ops if op in "MIX=S")
        start, end = left, left + query
        if flag & 0x10:
            start, end = lengths[r[0]] - end, lengths[r[0]] - start
        out[r[0]].append((start, end, r[2], r[5]))
    return out


def split_intervals(sam, piece):
    out = collections.defaultdict(list)
    for line in open(sam):
        if line.startswith("@"):
            continue
        r = line.split("\t", 10)
        flag, mapq = int(r[1]), int(r[4])
        if flag & 0x100 or mapq < 4:
            continue
        m = re.match(r"(.*)_p(\d+)$", r[0])
        i = int(m.group(2))
        out[m.group(1)].append((i * piece, (i + 1) * piece, r[2], r[5]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--piece", type=int, default=1000)
    args = ap.parse_args()
    w = args.work
    truth = read_table(f"{w}/hifi.truth.tsv")
    lengths = {}
    with open(f"{w}/hifi.fq") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].strip()
            elif i % 4 == 1:
                lengths[name] = len(line.strip())
    markers = collections.defaultdict(list)
    for row in read_table(args.markers):
        markers[(row["accession"], row["contig"])].append((int(row["start"]), int(row["end"]), row["marker"]))
    native = native_intervals(f"{w}/out/native.sam", lengths)
    split = split_intervals(f"{w}/out/split{args.piece}.sam", args.piece)

    counts = collections.Counter()
    examples = []
    for row in truth:
        name, start, end = row["read"], int(row["start"]), int(row["end"])
        for m_start, m_end, marker in markers[(row["accession"], row["contig"])]:
            overlap = min(end, m_end) - max(start, m_start)
            if overlap < 200:
                continue
            # the marker on the read, forward read coordinates (reads carry indels: approximate)
            a, b = max(m_start, start) - start, min(m_end, end) - start
            if row["strand"] == "-":
                a, b = (end - start) - b, (end - start) - a
            cut = m_start < start or m_end > end
            in_native = any(min(b, e) - max(a, s) >= 0.5 * (b - a) for s, e, _, _ in native.get(name, []))
            in_split = any(min(b, e) - max(a, s) > 0 for s, e, _, _ in split.get(name, []))
            counts[("found" if in_native else "missed", "cut" if cut else "whole", "split finds it" if in_split else "split misses it")] += 1
            if not in_native and in_split and len(examples) < 8:
                examples.append((name, row["strand"], lengths[name], marker, a, b, cut,
                                 sorted(native.get(name, []))))
    for key, n in sorted(counts.items()):
        print("\t".join(key), n, sep="\t")
    for e in examples:
        print("\nEXAMPLE", *e[:7])
        for rec in e[7]:
            print("   native record", rec)


if __name__ == "__main__":
    main()
