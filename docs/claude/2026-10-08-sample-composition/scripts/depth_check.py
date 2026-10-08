#!/usr/bin/env python3
"""Where the composition's error comes from: per called species, protal's depth (.profile.log VCov) against the depth its
simulated reads give (the manifest's vertical coverage, read bases over genome length, summed over the species' strains,
times the sample's FragmentBaseShare: a pair's overlap counted once, as in protal's depth), and the genome size protal
takes (GenomeSize, the species' mean in species_priors.tsv) against the length of the genomes simulated (weighted by
their coverage). Prints per run the median and the 10-90% range of both ratios over the species, and the ratio of the
bases explained (depth x size summed) to the bases of those species' reads.

Usage: depth_check.py --manifest sims/manifest.tsv --run all=run_db_all [--run held=run_db_held]
"""
import argparse
import collections
import csv
import glob
import os
import statistics


def quantiles(values):
    values = sorted(values)
    pick = lambda q: values[min(len(values) - 1, int(q * len(values)))]
    return f"median {statistics.median(values):.4f} (10-90%: {pick(0.1):.4f}-{pick(0.9):.4f}, n={len(values)})"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--run", action="append", required=True, help="NAME=OUTPUT_DIR of a protal run")
    args = ap.parse_args()
    truth = collections.defaultdict(lambda: collections.defaultdict(lambda: [0.0, 0.0, 0]))  # sample -> species -> [cov, cov*len, pairs]
    with open(args.manifest) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            t = truth[row["sample"]][row["taxonomy"].split(";")[-1]]
            cov = float(row["vertical_coverage"])
            t[0] += cov
            t[1] += cov * int(row["genome_length"])
            t[2] += int(row["read_pairs"])
    for spec in args.run:
        name, folder = spec.split("=", 1)
        depth_ratios, size_ratios, explained = [], [], []
        for sample, species in sorted(truth.items()):
            log = glob.glob(os.path.join(folder, "**", f"{sample}.profile.log"), recursive=True)[0]
            with open(log[: -len(".log")] + ".composition") as fh:
                share = float(next(csv.DictReader(fh, delimiter="\t"))["FragmentBaseShare"])
            protal_bases = true_bases = 0.0
            with open(log) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    if row["Predicted"] != "1" or row["Name"] not in species:
                        continue
                    cov, cov_len, _ = species[row["Name"]]
                    depth, size = float(row["VCov"]), float(row["GenomeSize"])
                    depth_ratios.append(depth / (cov * share))
                    size_ratios.append(size / (cov_len / cov))
                    protal_bases += depth * size
                    true_bases += cov_len * share
            explained.append(protal_bases / true_bases)
        print(f"{name}: depth / simulated depth {quantiles(depth_ratios)}")
        print(f"{name}: genome size / simulated genomes' {quantiles(size_ratios)}")
        print(f"{name}: bases explained / the called species' read bases, per sample {quantiles(explained)}")


if __name__ == "__main__":
    main()
