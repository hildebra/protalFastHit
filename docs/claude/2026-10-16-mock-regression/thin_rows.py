#!/usr/bin/env python3
"""How the simulated training rows of a build look at the evidence of the mocks' false calls: per bin of fragments and
top identity (>= 0.97 or below), the rows, the share present (truth 1) and the share the model called (prediction 1),
for the rows of all samples and of deep samples (meta_read_pairs >= DEEP). A row is a taxon with reads in a simulated
sample: every read comes from a member genome (its species, an in-silico strain, a held-out relative) or the host.

    python3 thin_rows.py TRAINING_DATA.tsv          prints the table
"""
import collections
import csv
import sys

BINS = [(1, 2), (3, 5), (6, 10), (11, 25), (26, 100), (101, 10 ** 12)]
DEEP = 5_000_000


def label(n):
    return next(f"{a}-{b}" if b < 10 ** 12 else f">{a - 1}" for a, b in BINS if a <= n <= b)


def main(path):
    csv.field_size_limit(1 << 30)
    counts = collections.Counter()
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            try:
                n = round(float(r["fragments"]))
                identity = float(r["top_identity"])
                pairs = float(r["meta_read_pairs"])
            except (ValueError, KeyError):
                continue
            if n < 1:
                continue
            high = "high" if identity >= 0.97 else "low"
            for depth in ("all", "deep") if pairs >= DEEP else ("all",):
                key = (depth, label(n), high)
                counts[key + ("rows",)] += 1
                counts[key + ("present",)] += r["truth"] in ("1", "1.0", "TRUE", "True")
                counts[key + ("called",)] += r["prediction"] in ("1", "1.0", "TRUE", "True")
    print("samples\tfragments\ttop_identity\trows\tpresent_share\tcalled_share")
    for depth in ("all", "deep"):
        for a, b in BINS:
            for high in ("high", "low"):
                key = (depth, f"{a}-{b}" if b < 10 ** 12 else f">{a - 1}", high)
                rows = counts[key + ("rows",)]
                if rows:
                    print(f"{depth}\t{key[1]}\t{'>=0.97' if high == 'high' else '<0.97'}\t{rows}\t"
                          f"{counts[key + ('present',)] / rows:.3f}\t{counts[key + ('called',)] / rows:.3f}")


if __name__ == "__main__":
    main(sys.argv[1])
