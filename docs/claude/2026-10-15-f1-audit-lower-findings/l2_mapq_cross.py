#!/usr/bin/env python3
"""The primary records of two runs of the same samples, mate by mate (F1 audit L2: --gene_end_clips off against
credit): how many mates changed MAPQ (up, down) or taxon, and how many crossed the profiler's cutoff (MAPQ 4), in
each direction, on the same taxon or to another one. A mate "counts" with MAPQ 4 or more.

  columns: mates (with a record in either run), taxon_changed, mapq_up, mapq_down (same taxon), counted_lost and
  counted_gained (same taxon: the MAPQ crossed 4), counted_moved (counted in both on two taxa), counted_lost_moved
  and counted_gained_moved (crossed 4 with the taxon changed).

usage: l2_mapq_cross.py OUT.tsv SAMPLE:SAM_A:SAM_B [...]
Standard library only (zstdcat for .zst SAMs).
"""
import collections
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "2026-10-15-mate-consensus"))
from compare_records import CUTOFF, read_records  # noqa: E402

COLUMNS = ["mates", "taxon_changed", "mapq_up", "mapq_down", "counted_lost", "counted_gained", "counted_moved",
           "counted_lost_moved", "counted_gained_moved"]


def main():
    out, specs = sys.argv[1], sys.argv[2:]
    rows = []
    for spec in specs:
        name, path_a, path_b = spec.split(":")
        ra, rb = read_records(path_a), read_records(path_b)
        c = collections.Counter()
        for qname in set(ra) | set(rb):
            ea = ra.get(qname, [None] * 4)
            eb = rb.get(qname, [None] * 4)
            for m in (0, 2):
                ta, qa, tb, qb = ea[m], ea[m + 1], eb[m], eb[m + 1]
                if ta is None and tb is None:
                    continue
                c["mates"] += 1
                ca = ta is not None and qa >= CUTOFF
                cb = tb is not None and qb >= CUTOFF
                if ta == tb:
                    c["mapq_up"] += qb > qa
                    c["mapq_down"] += qb < qa
                    c["counted_lost"] += ca and not cb
                    c["counted_gained"] += cb and not ca
                else:
                    c["taxon_changed"] += 1
                    c["counted_moved"] += ca and cb
                    c["counted_lost_moved"] += ca and not cb
                    c["counted_gained_moved"] += cb and not ca
        row = {"sample": name}
        row.update({k: c.get(k, 0) for k in COLUMNS})
        rows.append(row)
        print(name, " ".join(f"{k}={row[k]}" for k in COLUMNS), file=sys.stderr, flush=True)
    total = {"sample": "all"}
    total.update({k: sum(r[k] for r in rows) for k in COLUMNS})
    rows.append(total)
    with open(out, "w") as fh:
        fh.write("\t".join(["sample"] + COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[k]) for k in ["sample"] + COLUMNS) + "\n")


if __name__ == "__main__":
    main()
