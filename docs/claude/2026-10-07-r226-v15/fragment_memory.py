#!/usr/bin/env python3
"""The memory error_reads.py's extract() holds per tracked fragment, as r226 v15 ran it ("v15": pass 1's why, qname ->
set of reasons, and pass 2's fragments, qname -> a Fragment with every aligned and seeded-on taxon in sets) and as
2543408 does ("now": one dict, qname -> a tuple of shared reason strings, then its Fragment with only the sample's
error taxa), from synthetic paired-end records (ART-like names, two records a fragment, --zf taxids on each record's
ZF tag, --error-taxa of the 5,000 taxa the sample's errors). Prints the resident set's growth per fragment.

    python3 fragment_memory.py --fragments 1000000 --zf 8
"""
import argparse
import collections
import random

NONE = frozenset()


class V15Fragment:
    __slots__ = ("aligned", "best", "best_mapq", "failed", "records")

    def __init__(self):
        self.aligned, self.best, self.best_mapq, self.failed, self.records = set(), None, 0, set(), 0


class Fragment:  # error_reads.Fragment since 2543408
    __slots__ = ("reasons", "aligned", "best", "best_mapq", "failed", "records")

    def __init__(self, reasons=()):
        self.reasons, self.aligned, self.best, self.best_mapq, self.failed, self.records = reasons, NONE, None, 0, NONE, 0


def rss():
    with open("/proc/self/status") as fh:
        for line in fh:
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024
    return 0


def fields(line):
    f = line.split(b"\t", 5)
    zf = line.rfind(b"\tZF:Z:")
    return f[0], f[2], line[zf + 6:].rstrip(b"\n").split(b",") if zf >= 0 else ()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--fragments", type=int, default=1_000_000)
    ap.add_argument("--zf", type=int, default=8, help="taxids on a record's ZF tag")
    ap.add_argument("--error-taxa", type=int, default=50, help="of the 5,000 taxa, the sample's error taxa")
    ap.add_argument("--layout", choices=("v15", "now"), default="now")
    opts = ap.parse_args()
    rng = random.Random(1)
    taxa = [str(rng.randrange(1, 140000)).encode() for _ in range(5000)]
    error_taxa = set(taxa[:opts.error_taxa])
    lines = []
    for i in range(opts.fragments):
        qname = f"GCA_{rng.randrange(10**9):09d}.1_contig{rng.randrange(1000)}-{i}".encode()
        for _ in range(2):
            zf = b",".join(rng.choice(taxa) for _ in range(opts.zf))
            lines.append(qname + b"\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\tZU:i:0\tZT:i:0\tZF:Z:" + zf + b"\n")
    reason = "source:12345"
    base = rss()
    if opts.layout == "v15":
        why = collections.defaultdict(set)
        for line in lines:
            qname, _, _ = fields(line)
            why[qname].add(reason)
        fragments = {}
        for line in lines:
            qname, _, failed = fields(line)
            fragment = fragments.get(qname)
            if fragment is None:
                fragment = fragments[qname] = V15Fragment()
            fragment.records += 1
            fragment.failed.update(failed)
    else:
        why, names, shared = {}, {}, {}
        for line in lines:
            qname, _, _ = fields(line)
            r = names.setdefault(reason, reason)
            have = why.get(qname, ())
            if r not in have:
                why[qname] = have + (r,)
        for line in lines:
            qname, _, failed = fields(line)
            fragment = why[qname]
            if isinstance(fragment, tuple):
                fragment = why[qname] = Fragment(fragment)
            fragment.records += 1
            for taxid in failed:
                if taxid in error_taxa:
                    if fragment.failed is NONE:
                        fragment.failed = set()
                    fragment.failed.add(shared.setdefault(taxid, taxid))
    grown = rss() - base
    print(f"{opts.layout}: {opts.fragments} fragments, {opts.zf} ZF taxa a record, {opts.error_taxa} error taxa: "
          f"{grown / opts.fragments:.0f} B per fragment ({grown / 1e9:.2f} GB)")


if __name__ == "__main__":
    main()
