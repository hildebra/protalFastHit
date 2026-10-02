#!/usr/bin/env python3
"""Where gene_neighbours.py placed the database genes in each genome (gene_positions.tsv) against where the
simulator put them (simulate_gtdb_release.py's simulation/marker_positions.tsv): per genome kind
(representative or other strain) and placement (exact or trace), how many are at the true place (same contig and
strand, start within --tolerance), elsewhere, or missing.

usage: check_placements.py --db CONVERTED_DB --world WORLD [--tolerance 30]
"""

import argparse
import collections
import csv
import os


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--db", required=True)
    ap.add_argument("--world", required=True)
    ap.add_argument("--tolerance", type=int, default=30)
    opts = ap.parse_args()
    with open(os.path.join(opts.db, "gene2geneid.tsv")) as fh:
        gene_id = {f[0]: int(f[1]) for f in (line.rstrip("\n").split("\t") for line in fh) if len(f) >= 2}
    reps = set()
    with open(os.path.join(opts.db, "internal_taxonomy.dmp")) as fh:
        next(fh)
        for f in (line.rstrip("\n").split("\t") for line in fh):
            if f[4] == "species" and len(f) > 6:
                reps.add(f[6])
    truth = {}
    with open(os.path.join(opts.world, "simulation", "marker_positions.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["marker"] in gene_id:
                truth.setdefault((r["accession"], gene_id[r["marker"]]), (r["contig"], int(r["start"]), int(r["end"]), r["strand"]))
    counts = collections.Counter()
    placed = set()
    with open(os.path.join(opts.db, "gene_positions.tsv")) as fh:
        rows = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        for r in rows:
            key = (r["accession"], int(r["gene"]))
            placed.add(key)
            kind = "representative" if r["accession"] in reps else "strain"
            t = truth.get(key)
            if t is None:
                verdict = "not in the simulator's genome"
            elif t[0] == r["contig"] and t[3] == r["strand"] and abs(t[1] - int(r["start"])) <= opts.tolerance:
                verdict = "at its place" if t[1] == int(r["start"]) and t[2] == int(r["end"]) else "near its place"
            else:
                verdict = "elsewhere"
            counts[(kind, r["placed"], verdict)] += 1
    genomes = {a for a, _ in placed}
    for (acc, gene) in truth:
        if acc in genomes and (acc, gene) not in placed:
            counts[("representative" if acc in reps else "strain", "-", "missed")] += 1
    for key in sorted(counts):
        print("\t".join(key), counts[key], sep="\t")


if __name__ == "__main__":
    main()
