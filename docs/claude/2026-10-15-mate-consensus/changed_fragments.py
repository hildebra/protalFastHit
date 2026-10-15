#!/usr/bin/env python3
"""The simulated fragments whose counted species changed between two runs (compare_records.py's "moved" fragments),
each with its truth and both runs' primary records (species, gene, position, CIGAR, MAPQ), for a look at what moved
them; and a tally by the kind of change.

  kinds: right_lost  (the run A counted the fragment for its own species, B does not: B's MAPQ under 4, or another
                      species)
         wrong_gained (B counts it for a species A did not, not its own)
         right_gained, wrong_lost: the converse

usage: changed_fragments.py OUT.tsv DATABASE_FOLDER WORK SAMPLE:SAM_A:SAM_B [...]
Standard library only (zstdcat for .zst SAMs).
"""
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_records import CUTOFF, load_db_species, load_truth, open_sam  # noqa: E402


def read_full(path):
    """qname -> {mate: (taxid, gene, pos, cigar, mapq)} of the primary records."""
    fragments = collections.defaultdict(dict)
    for line in open_sam(path):
        if line.startswith("@"):
            continue
        f = line.split("\t", 7)
        flag = int(f[1])
        if flag & 0x904:
            continue
        taxid, gene = f[2].split("_")[:2]
        fragments[f[0]][2 if flag & 0x80 else 1] = (taxid, gene, f[3], f[5], int(f[4]))
    return fragments


def counted(records, species):
    return {species.get(r[0]) for r in records.values() if r[4] >= CUTOFF}


def main():
    out, database, work = sys.argv[1:4]
    species = load_db_species(database)
    truth = load_truth(work)
    tally = collections.Counter()
    with open(out, "w") as fh:
        fh.write("sample\tqname\tkind\ttruth\trecords_a\trecords_b\n")
        for spec in sys.argv[4:]:
            name, path_a, path_b = spec.split(":")
            ra, rb = read_full(path_a), read_full(path_b)
            for qname in set(ra) | set(rb):
                a, b = ra.get(qname, {}), rb.get(qname, {})
                if a == b:
                    continue
                own = truth.get(qname.rsplit("-", 1)[0])
                ca, cb = counted(a, species), counted(b, species)
                kinds = []
                if own in ca and own not in cb:
                    kinds.append("right_lost")
                if own in cb and own not in ca:
                    kinds.append("right_gained")
                if cb - ca - {own}:
                    kinds.append("wrong_gained")
                if ca - cb - {own}:
                    kinds.append("wrong_lost")
                for k in kinds:
                    tally[k] += 1
                if not kinds:
                    continue
                show = lambda recs: "; ".join(f"m{m}:{species.get(r[0], r[0]).replace('s__', '')}|{r[0]}_{r[1]}@{r[2]} "
                                              f"{r[3]} q{r[4]}" for m, r in sorted(recs.items())) or "-"
                fh.write(f"{name}\t{qname}\t{','.join(kinds)}\t{(own or '?').replace('s__', '')}\t{show(a)}\t{show(b)}\n")
    for k in ("right_lost", "right_gained", "wrong_gained", "wrong_lost"):
        print(f"{k}\t{tally[k]}")


if __name__ == "__main__":
    main()
