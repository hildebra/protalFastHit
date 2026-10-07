#!/usr/bin/env python3
"""simulate_metagenomes --long_samples on a community of real-size genomes: wall and CPU time of whole samples
(template drawing, read model, compression) per read setup and thread count, and the cost of a genome read (a sample
of many genomes and few bases), for the collector's estimates (LONG_SECONDS_PER_MB, LONG_GENOME_SECONDS).

    bench_community.py --simulator BIN --model QSHMM-ONT-HQ.model --genomes DIR --out DIR [--threads 1,4]
"""
import argparse
import glob
import os
import random
import resource
import subprocess
import sys
import time


def run(command):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    began = time.perf_counter()
    done = subprocess.run(command, capture_output=True, text=True)
    wall = time.perf_counter() - began
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    if done.returncode != 0:
        sys.exit(f"{' '.join(command)} failed: {done.stderr[-2000:]}")
    return wall, after.ru_utime - before.ru_utime + after.ru_stime - before.ru_stime, after.ru_maxrss


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--simulator", required=True)
    p.add_argument("--model", required=True)
    p.add_argument("--genomes", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--threads", default="1,4")
    a = p.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fastas = sorted(glob.glob(os.path.join(a.genomes, "*.fna.gz")))
    rng = random.Random(3)
    weights = [rng.lognormvariate(0, 1.5) * 3.4e6 for _ in fastas]  # lognormal abundances, sigma 1.5

    def tasks(name, bases, genomes):
        samples = os.path.join(a.out, name + "_samples.tsv")
        table = os.path.join(a.out, name + "_genomes.tsv")
        with open(samples, "w") as fh:
            fh.write(f"sample\tout\tbases\tseed\n{name}\t{os.path.join(a.out, name + '.fq.zst')}\t{bases}\t7\n")
        with open(table, "w") as fh:
            fh.write("sample\tgenome\tfasta\tweight\thost\n")
            for i, (fasta, w) in enumerate(zip(fastas[:genomes], weights)):
                fh.write(f"{name}\tG{i}\t{fasta}\t{w!r}\t0\n")
        return samples, table

    setups = [("hifi", "hifi:15000:3000:3", 1_000_000_000), ("ultima", "ultima:300:40:25:2", 500_000_000),
              ("ont", "qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36", 1_000_000_000)]
    print("run\tsetup\tgenomes\tGb\tthreads\twall s\tCPU s\tCPU s per Mb\tmax RSS MB")
    for kind, setup, bases in setups:
        samples, table = tasks(kind, bases, len(fastas))
        for threads in a.threads.split(","):
            stats = os.path.join(a.out, f"{kind}_{threads}.stats")
            wall, cpu, rss = run([a.simulator, "--long_samples", samples, "--long_genomes", table, "--long_setup", setup,
                                  "--long_model", a.model, "--long_stats", stats, "-t", threads])
            print(f"{kind}\t{setup}\t{len(fastas)}\t{bases / 1e9:g}\t{threads}\t{wall:.1f}\t{cpu:.1f}\t"
                  f"{cpu / (bases / 1e6):.4f}\t{rss / 1024:.0f}", flush=True)
    # The genome reads: every genome gets about one 1 kb read (equal weights, as many reads as genomes).
    samples, table = tasks("genomes", 1000 * len(fastas), len(fastas))
    with open(table) as fh:
        lines = fh.readlines()
    with open(table, "w") as fh:
        fh.write(lines[0] + "".join("\t".join(line.split("\t")[:3] + ["1.0", "0\n"]) for line in lines[1:]))
    stats = os.path.join(a.out, "genomes.stats")
    wall, cpu, rss = run([a.simulator, "--long_samples", samples, "--long_genomes", table, "--long_setup",
                          "hifi:1000:0:3", "--long_stats", stats, "-t", "1"])
    print(f"genome reads\thifi:1000:0:3\t{len(fastas)}\t{1000 * len(fastas) / 1e9:g}\t1\t{wall:.1f}\t{cpu:.1f}\t"
          f"{cpu / len(fastas):.4f} s per genome (~{1 - (1 - 1 / len(fastas)) ** len(fastas):.2f} of them read)\t"
          f"{rss / 1024:.0f}")


if __name__ == "__main__":
    main()
