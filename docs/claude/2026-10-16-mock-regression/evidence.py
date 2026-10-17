#!/usr/bin/env python3
"""The evidence behind each call of the mock-community runs in a score_mocks.py results folder (an unpacked
--store_results_tar): per run, sample and species called or expected, its status (scores/species.tsv), and from the
sample's .profile.log (the Summary column; 0.6.0a's file has no header, same column order) its probability, genes,
fragments (Total Hits), mean ANI, mean MAPQ and depth. A false species is marked by whether its genus has a member in
the sample's truth (a relative's reads) and whether it is a member of another truth table (another mock's species).

    python3 evidence.py RESULTS_DIR OUT_DIR          writes OUT_DIR/calls.tsv
"""
import csv
import glob
import os
import re
import sys

SUMMARY = re.compile(r"Genes: (\d+), HittableGenes: (\d+), Total Hits: (\d+), Total MGenome: (\d+), "
                     r"MeanANI: ([-\d.]+), MeanMAPQ: ([-\d.]+), VCOV: ([-\d.]+)")
COLUMNS = ["run", "set", "sample", "read_type", "truth", "species", "status", "truth_share", "profile_share",
           "probability", "genes", "hittable_genes", "fragments", "mean_ani", "mean_mapq", "vcov", "genus_in_truth",
           "in_other_truth"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def genus(species):
    return species.removeprefix("s__").split(" ")[0]


def evidence(path):
    """{species: (probability, genes, hittable, fragments, ani, mapq, vcov)} of a .profile.log."""
    found = {}
    with open(path, newline="") as fh:
        for f in csv.reader(fh, delimiter="\t"):
            if len(f) < 11 or f[0] == "Predicted":
                continue
            m = SUMMARY.search(f[10])
            if not m:
                continue
            name = f[3].split(";")[-1]
            found[name] = (f[1],) + m.group(1, 2, 3, 5, 6, 7)
    return found


def main(results, out):
    truths = {}
    for path in glob.glob(os.path.join(results, "inputs", "truth", "*.tsv")):
        rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip() and not l.startswith("#")][1:]
        truths[os.path.basename(path)[:-4]] = {r[0] for r in rows if r[0] != "?" and r[1] != "0"}
    samples = {r["sample"]: r for r in read_tsv(os.path.join(results, "inputs", "samples.tsv"))}
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "calls.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for run in sorted(os.listdir(results)):
            species_path = os.path.join(results, run, "scores", "species.tsv")
            if not os.path.isfile(species_path):
                continue
            logs = {os.path.basename(p)[:-len(".profile.log")]: p
                    for p in glob.glob(os.path.join(results, run, "**", "*.profile.log"), recursive=True)}
            cache = {}
            for r in read_tsv(species_path):
                s = samples[r["sample"]]
                if r["sample"] not in cache:
                    cache[r["sample"]] = evidence(logs[r["sample"]]) if r["sample"] in logs else {}
                e = cache[r["sample"]].get(r["species"], ("NA",) * 7)
                own = truths.get(s["truth"], set())
                others = set().union(*(v for k, v in truths.items() if k != s["truth"])) - own
                w.writerow([run, r["set"], r["sample"], s["read_type"], s["truth"], r["species"], r["status"],
                            r["truth_share"], r["profile_share"], *e,
                            int(genus(r["species"]) in {genus(k) for k in own}), int(r["species"] in others)])


if __name__ == "__main__":
    main(*sys.argv[1:3])
