#!/usr/bin/env python3
"""Would a per-gene correction of the gene-end loss fix the depth? On run_validation.sh's db_all run (every species in the
database), from each called species' genes in <sample>.profile.genes.log (VCov: fragment bases over gene length, as the
depth's per-gene values): c, the bases a gene loses at its ends, fitted as the median of (1 - gene depth / true depth) x
length over genes with 20 or more reads; then each species' depth as the median of its genes' depths, uncorrected and
divided by (1 - c / length), against the depth its reads give (manifest coverage x FragmentBaseShare), per domain; and
the share of each sample's bases the species explain with either depth. Run in the validation folder.

Usage: edge_correction.py [--run run_db_all] [--manifest sims/manifest.tsv]
"""
import argparse
import collections
import csv
import glob
import statistics


def spread(values):
    values = sorted(values)
    pick = lambda q: values[min(len(values) - 1, int(q * len(values)))]
    return f"median {statistics.median(values):.4f} (10-90%: {pick(0.1):.4f}-{pick(0.9):.4f}, n={len(values)})"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--run", default="run_db_all")
    ap.add_argument("--manifest", default="sims/manifest.tsv")
    args = ap.parse_args()
    truth = collections.defaultdict(lambda: collections.defaultdict(lambda: [0.0, 0.0]))  # coverage, coverage x length
    with open(args.manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            t = truth[row["sample"]][row["taxonomy"].split(";")[-1]]
            t[0] += float(row["vertical_coverage"])
            t[1] += float(row["vertical_coverage"]) * int(row["genome_length"])
    genes = []  # (sample, species, domain, length, gene depth / true depth, reads)
    composition = {}
    for sample in truth:
        comp = glob.glob(f"{args.run}/**/{sample}.profile.composition", recursive=True)[0]
        composition[sample] = next(csv.DictReader(open(comp), delimiter="\t"))
        share = float(composition[sample]["FragmentBaseShare"])
        log = glob.glob(f"{args.run}/**/{sample}.profile.genes.log", recursive=True)[0]
        with open(log) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                name = row["Lineage"].split(";")[-1]
                if row["Predicted"] != "1" or name not in truth[sample]:
                    continue
                true_depth = truth[sample][name][0] * share
                genes.append((sample, name, row["Lineage"].split(";")[0], int(row["GeneRefLength"]),
                              float(row["VCov"]) / true_depth, int(row["TotalReads"])))
    deep = [(1 - ratio) * length for _, _, _, length, ratio, reads in genes if reads >= 20]
    c = statistics.median(deep)
    print(f"c (bases lost per gene): {c:.1f}, from {len(deep)} genes with 20 or more reads")
    for label, b in (("by gene length, 0-299", (0, 300)), ("300-599", (300, 600)), ("600-1199", (600, 1200)),
                     ("1200 or more", (1200, 10 ** 9))):
        raw = [r for _, _, _, l, r, n in genes if b[0] <= l < b[1] and n >= 20]
        fixed = [r / (1 - c / l) for _, _, _, l, r, n in genes if b[0] <= l < b[1] and n >= 20]
        if raw:
            print(f"  genes {label}: depth ratio {statistics.median(raw):.4f} -> corrected {statistics.median(fixed):.4f}")
    by_species = collections.defaultdict(list)
    for sample, name, domain, length, ratio, _ in genes:
        by_species[(sample, name, domain)].append((ratio, ratio / (1 - c / length)))
    for domain in sorted({k[2] for k in by_species}):
        raw = [statistics.median(r for r, _ in v) for k, v in by_species.items() if k[2] == domain]
        fixed = [statistics.median(f for _, f in v) for k, v in by_species.items() if k[2] == domain]
        print(f"{domain}: species depth / true, uncorrected {spread(raw)}")
        print(f"{domain}: species depth / true, corrected   {spread(fixed)}")
    for sample in sorted(truth):
        comp = composition[sample]
        total = float(comp["ScannedBases"]) * float(comp["FragmentBaseShare"])
        raw = fixed = 0.0
        for (s, name, _), v in by_species.items():
            if s != sample:
                continue
            length = truth[s][name][1] / truth[s][name][0]
            true_depth = truth[s][name][0] * float(comp["FragmentBaseShare"])
            raw += statistics.median(r for r, _ in v) * true_depth * length
            fixed += statistics.median(f for _, f in v) * true_depth * length
        print(f"{sample}: explained by the median gene depths {raw / total:.4f}, corrected {fixed / total:.4f} "
              f"(protal: {float(comp['AttributedShare']):.4f})")


if __name__ == "__main__":
    main()
