#!/usr/bin/env python3
"""Does the collector's thread pool run long-read samples side by side? long_read_sample's Python part for PacBio
HiFi (long_read_templates, then hifi_reads.simulate) for N samples of one community (a sim_cost.sh manifest):
one after the other, on N threads (concurrent.futures.ThreadPoolExecutor, as simulate_long does), and in N
processes (ProcessPoolExecutor).

usage: gil_check.py MANIFEST OUT [--samples 4] [--bases 30000000]
"""
import argparse
import concurrent.futures
import csv
import os
import sys
import time

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import collect_training_data as ctd  # noqa: E402
import hifi_reads  # noqa: E402


def one(task):
    names, error = ctd.long_read_templates(task)
    if error:
        raise RuntimeError(error)
    hifi_reads.simulate(task["templates"], task["out"], task["setup"]["q_sd"], task["seed"])
    return len(names)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest")
    p.add_argument("out")
    p.add_argument("--samples", type=int, default=4)
    p.add_argument("--bases", type=int, default=30_000_000)
    opts = p.parse_args()
    os.makedirs(opts.out, exist_ok=True)
    with open(opts.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    first = rows[0]["sample"]
    genomes = [{"genome": r["genome"], "fasta": r["fasta_path"],
                "weight": float(r["relative_abundance"]) * float(r["genome_length"])} for r in rows if r["sample"] == first]
    setup = ctd.parse_long_setup("hifi:15000:3000:3")
    tasks = [{"seed": 100 + i, "genomes": genomes, "setup": setup, "bases": opts.bases,
              "templates": os.path.join(opts.out, f"t{i}.fa"), "out": os.path.join(opts.out, f"r{i}.fq.gz")}
             for i in range(opts.samples)]
    print(f"mode\tsamples\twall_s\t{len(genomes)} genomes, {opts.bases} bases each", flush=True)
    began = time.time()
    for task in tasks:
        one(task)
    print(f"serial\t{opts.samples}\t{time.time() - began:.1f}", flush=True)
    for name, pool in (("threads", concurrent.futures.ThreadPoolExecutor),
                       ("processes", concurrent.futures.ProcessPoolExecutor)):
        began = time.time()
        with pool(opts.samples) as executor:
            list(executor.map(one, tasks))
        print(f"{name}\t{opts.samples}\t{time.time() - began:.1f}", flush=True)


if __name__ == "__main__":
    main()
