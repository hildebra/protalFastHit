#!/usr/bin/env python3
"""Compares protal's sample compositions (<profile>.composition, <profile>) with a simulate_metagenomes truth.

Per sample and run (a database): the share of the read pairs that come from species the database holds, and from the
species protal called correctly (what the attributed share should be given the calls); the truth's unknown share, the
cells of species the database lacks (the simulator's relative abundances are cell shares: coverage, normalised), and
of species it holds that protal did not call; the average genome size of all cells and of the cells of the correctly
called species; the species the database lacks. Against protal's AttributedShare, UnknownShare (the profile's "?"),
AverageGenomeSize and MissingSpeciesAt*Depth.

Usage: evaluate.py --manifest sims/manifest.tsv --held_out held_out.txt --run all=run_db_all --run held=run_db_held
"""
import argparse
import collections
import csv
import glob
import os
import sys


def read_manifest(path):
    """{sample: [(species, genome_length, read_pairs, relative_abundance)]}."""
    samples = collections.defaultdict(list)
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            species = row["taxonomy"].split(";")[-1]
            samples[row["sample"]].append((species, int(row["genome_length"]), int(row["read_pairs"]),
                                           float(row["relative_abundance"])))
    return samples


def read_profile(path):
    """({called species name}, unknown share or None) of a .profile."""
    called, unknown = set(), None
    with open(path) as fh:
        for line in fh:
            rep, lineage, share = line.rstrip("\n").split("\t")
            if lineage == "?":
                unknown = float(share)
            else:
                called.add(lineage.split(";")[-1])
    return called, unknown


def read_composition(path):
    with open(path) as fh:
        return next(csv.DictReader(fh, delimiter="\t"))


def number(text):
    return float("nan") if text in ("NA", "") else float(text)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--held_out", required=True, help="species the held-out database lacks, one per line")
    ap.add_argument("--run", action="append", required=True, help="NAME=OUTPUT_DIR of a protal run; NAME 'held' lacks --held_out")
    args = ap.parse_args()
    truth = read_manifest(args.manifest)
    with open(args.held_out) as fh:
        held_out = {line.strip() for line in fh if line.strip()}

    out = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    out.writerow(["run", "sample", "species", "lacking", "called", "true_called", "pairs_in_db", "pairs_true_called",
                  "attributed_share", "unknown_true_db", "unknown_true_calls", "unknown_profile", "ags_all", "ags_true_called",
                  "ags_protal", "ags_p10", "ags_p90", "missing_at_median", "missing_at_lowest"])
    for spec in args.run:
        name, folder = spec.split("=", 1)
        lacking = held_out if name == "held" else set()
        for sample, genomes in sorted(truth.items()):
            profiles = glob.glob(os.path.join(folder, "**", f"{sample}.profile"), recursive=True)
            if len(profiles) != 1:
                sys.exit(f"{folder}: {len(profiles)} profiles of {sample}")
            called, unknown = read_profile(profiles[0])
            comp = read_composition(profiles[0] + ".composition")
            present = {s for s, _, _, _ in genomes}
            true_called = called & (present - lacking)
            pairs = sum(p for _, _, p, _ in genomes)
            in_db = sum(p for s, _, p, _ in genomes if s not in lacking) / pairs
            on_calls = sum(p for s, _, p, _ in genomes if s in true_called) / pairs
            cells = sum(a for _, _, _, a in genomes)
            unknown_db = sum(a for s, _, _, a in genomes if s in lacking) / cells
            unknown_calls = sum(a for s, _, _, a in genomes if s not in true_called) / cells
            ags_all = sum(a * l for _, l, _, a in genomes) / cells
            called_cells = sum(a for s, _, _, a in genomes if s in true_called)
            ags_called = sum(a * l for s, l, _, a in genomes if s in true_called) / called_cells if called_cells else float("nan")
            out.writerow([name, sample, len(present), len(present & lacking), len(called), len(true_called),
                          f"{in_db:.4f}", f"{on_calls:.4f}", f"{number(comp['AttributedShare']):.4f}", f"{unknown_db:.4f}",
                          f"{unknown_calls:.4f}", "NA" if unknown is None else f"{unknown:.4f}", f"{ags_all:.0f}",
                          f"{ags_called:.0f}", comp["AverageGenomeSize"], comp["GenomeSizeP10"], comp["GenomeSizeP90"],
                          comp["MissingSpeciesAtMedianDepth"], comp["MissingSpeciesAtLowestDepth"]])


if __name__ == "__main__":
    main()
