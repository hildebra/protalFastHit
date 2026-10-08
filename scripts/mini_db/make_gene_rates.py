#!/usr/bin/env python3
"""make_gene_rates.py - the real genes' evolution speeds for simulate_gtdb_release.py --gene_rates r226.

Reads gene_congeners.tsv of a GTDB build (protal --build writes it into the database folder; build_gtdb_database.py
moves it to model_logs/): per gene, within_factor (how fast its copies diverge within species, over the species'
median gene: protal's conservation factor, GeneConservation.h) and between_factor (the same between congeneric
species). The gene ids are mapped to GTDB marker ids through the build's gene2geneid.tsv, or without one through
the converter's order (gtdb_to_protal_db.py: the bac120 markers sorted, then the ar53 markers not in bac120),
derived from the marker table.

Writes one line per marker and set (a marker in both sets has a line in each, with the same factors):
  set  marker  name  within_factor  between_factor  species
simulate_gtdb_release.py normalises each set's factors to a mean of 1 itself.

Usage:
  make_gene_rates.py --gene-congeners model_logs/gene_congeners.tsv [--gene-table gene2geneid.tsv]
      [--markers markers_r226.tsv] -o gene_rates_r226.tsv
"""

import argparse
import csv
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def read_markers(path):
    """[(set, marker, name)] of the marker table, in its order."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "set\t")) or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rows.append((f[0], f[1], f[2]))
    return rows


def converter_ids(markers):
    """{marker: gene id} as gtdb_to_protal_db.py numbers a release of these markers: the marker files sorted by set
    (bac120, then ar53) and by marker id, each marker numbered when first seen."""
    ids = {}
    for mset in ("bac120", "ar53"):
        for marker in sorted(m for s, m, _ in markers if s == mset):
            ids.setdefault(marker, len(ids) + 1)
    return ids


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--gene-congeners", required=True, help="gene_congeners.tsv of a GTDB build")
    ap.add_argument("--gene-table", help="that build's gene2geneid.tsv (default: the converter's order of --markers)")
    ap.add_argument("--markers", default=os.path.join(SCRIPT_DIR, "markers_r226.tsv"))
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args(argv)

    markers = read_markers(args.markers)
    if args.gene_table:
        with open(args.gene_table) as fh:
            ids = {f[0]: int(f[1]) for f in (line.split() for line in fh) if len(f) >= 2 and f[1].isdigit()}
    else:
        ids = converter_ids(markers)
    with open(args.gene_congeners) as fh:
        factors = {int(r["geneid"]): r for r in csv.DictReader(fh, delimiter="\t")}
    missing = sorted({m for _, m, _ in markers if m not in ids or ids[m] not in factors})
    if missing:
        sys.exit(f"{len(missing)} markers have no factors in {args.gene_congeners} ({', '.join(missing[:5])}): "
                 "another release's marker table, or a build of a gene subset?")
    if len(factors) != len(set(ids[m] for _, m, _ in markers)):
        sys.exit(f"{args.gene_congeners} has {len(factors)} genes, the marker table {len(set(ids.values()))}: "
                 "another release's genes?")
    with open(args.output, "w", newline="\n") as fh:
        fh.write(f"# per-gene evolution speeds of GTDB's markers, from {os.path.basename(args.gene_congeners)} "
                 "(make_gene_rates.py)\n")
        fh.write("set\tmarker\tname\twithin_factor\tbetween_factor\tspecies\n")
        for mset, marker, name in markers:
            r = factors[ids[marker]]
            fh.write(f"{mset}\t{marker}\t{name}\t{float(r['within_factor']):.4f}\t{float(r['between_factor']):.4f}\t"
                     f"{r['species']}\n")
    print(f"{len(markers)} markers of {len(factors)} genes: {args.output}")


if __name__ == "__main__":
    main()
