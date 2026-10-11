#!/usr/bin/env python3
"""Compares two protal SAMs (.sam.zst) record by record, ignoring @PG/@CO header lines: whether the record sets are
identical, and if not, which fields differ among the records of the reads whose records differ (per read, records
matched by reference and segment, else by order)."""
import argparse
import collections
import subprocess


def records(path):
    """{read name: [fields, ...]} of the alignment records, and the sorted header lines but @PG/@CO."""
    text = subprocess.run(["zstd", "-dcq", path], check=True, stdout=subprocess.PIPE).stdout.decode()
    by_read, header = collections.defaultdict(list), []
    for line in text.splitlines():
        if line.startswith("@"):
            if not line.startswith(("@PG", "@CO")):
                header.append(line)
            continue
        fields = line.split("\t")
        by_read[fields[0]].append(fields)
    return by_read, sorted(header)


def tags(fields):
    return {f[:2]: f[5:] for f in fields[11:]}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("a")
    p.add_argument("b")
    p.add_argument("--examples", type=int, default=3)
    args = p.parse_args()
    a, ha = records(args.a)
    b, hb = records(args.b)
    n_a = sum(len(v) for v in a.values())
    n_b = sum(len(v) for v in b.values())
    differing = [q for q in set(a) | set(b) if sorted(map(tuple, a.get(q, []))) != sorted(map(tuple, b.get(q, [])))]
    print(f"records: {n_a} vs {n_b}; reads: {len(a)} vs {len(b)}; header {'identical' if ha == hb else 'DIFFERS'}")
    if not differing:
        print("SAM records identical")
        return
    fields = collections.Counter()
    count_change = 0
    examples = []
    for q in differing:
        ra, rb = a.get(q, []), b.get(q, [])
        if len(ra) != len(rb):
            count_change += 1
            continue
        key = lambda f: (int(f[1]) & 0xC0, f[2], f[3])
        ra, rb = sorted(ra, key=key), sorted(rb, key=key)
        if [key(x) for x in ra] != [key(x) for x in rb]:
            ra, rb = a[q], b[q]
        what = set()
        for x, y in zip(ra, rb):
            for i, name in ((1, "FLAG"), (2, "RNAME"), (3, "POS"), (4, "MAPQ"), (5, "CIGAR"), (9, "SEQ")):
                if x[i] != y[i]:
                    what.add(name)
            tx, ty = tags(x), tags(y)
            for t in set(tx) | set(ty):
                if tx.get(t) != ty.get(t):
                    what.add("tag " + t)
        fields[" + ".join(sorted(what))] += 1
        if len(examples) < args.examples:
            examples.append((q, ra, rb))
    print(f"reads whose records differ: {len(differing)} of {len(set(a) | set(b))} "
          f"({100 * len(differing) / max(1, len(set(a) | set(b))):.3f}%); with another number of records: {count_change}")
    for what, n in fields.most_common():
        print(f"  {n:8d}  {what}")
    for q, ra, rb in examples:
        print(f"example {q}:")
        for x, y in zip(ra, rb):
            print("  a: " + "\t".join(x[1:6] + [f for f in x[11:] if not f.startswith(('MD', 'ZA'))])[:300])
            print("  b: " + "\t".join(y[1:6] + [f for f in y[11:] if not f.startswith(('MD', 'ZA'))])[:300])


if __name__ == "__main__":
    main()
