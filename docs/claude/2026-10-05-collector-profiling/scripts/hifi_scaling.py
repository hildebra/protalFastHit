#!/usr/bin/env python3
"""Does hifi_reads.simulate (numpy) run side by side on threads? The same N template files simulated one after the
other, on N threads of one process, and in N processes; Ultima (flow model) and HiFi reads.

usage: hifi_scaling.py --scripts PROTAL/scripts --genomes genomes.tsv --out DIR [--n 6]
"""

import argparse
import concurrent.futures
import os
import random
import resource
import sys
import time


def cpu(who=resource.RUSAGE_SELF):
    r = resource.getrusage(who)
    return r.ru_utime + r.ru_stime


def simulate(scripts, templates, out, setup):
    sys.path.insert(0, scripts)
    import hifi_reads
    return hifi_reads.simulate(templates, out, setup["q_sd"], 5, setup.get("q_mean"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--scripts", required=True)
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--bases", type=int, default=20_000_000)
    a = ap.parse_args()
    sys.path.insert(0, a.scripts)
    import collect_training_data as ctd
    os.makedirs(a.out, exist_ok=True)
    rows = [line.rstrip("\n").split("\t") for line in open(a.genomes)]
    rng = random.Random(3)
    genomes = [{"genome": n, "fasta": p, "weight": rng.lognormvariate(0, 1.5) * int(l)} for n, _, p, l in rows[:60]]
    for kind, text in (("ultima", "ultima:300:40:25:2"), ("pb", "hifi:15000:3000:3")):
        setup = ctd.parse_long_setup(text)
        files = []
        for i in range(a.n):
            path = os.path.join(a.out, f"{kind}{i}.fa")
            names, err = ctd.long_read_templates({"seed": 50 + i, "genomes": genomes, "setup": setup, "bases": a.bases,
                                                  "templates": path})
            assert not err
            files.append(path)
        outs = [os.path.join(a.out, f"{kind}{i}.fq.gz") for i in range(a.n)]
        results = {}
        w, c = time.time(), cpu()
        for f, o in zip(files, outs):
            simulate(a.scripts, f, o, setup)
        results["sequential"] = (time.time() - w, cpu() - c)
        w, c = time.time(), cpu()
        with concurrent.futures.ThreadPoolExecutor(a.n) as ex:
            list(ex.map(lambda fo: simulate(a.scripts, fo[0], fo[1], setup), zip(files, outs)))
        results["threads"] = (time.time() - w, cpu() - c)
        w, c = time.time(), cpu(resource.RUSAGE_CHILDREN)
        with concurrent.futures.ProcessPoolExecutor(a.n) as ex:
            list(ex.map(simulate, [a.scripts] * a.n, files, outs, [setup] * a.n))
        results["processes"] = (time.time() - w, cpu(resource.RUSAGE_CHILDREN) - c)
        base = results["sequential"][0]
        print(f"{kind}: {a.n} x {a.bases / 1e6:.0f} Mb: " + "; ".join(
            f"{mode} {wall:.1f} s wall, {c:.1f} s CPU ({base / wall:.2f}x)" for mode, (wall, c) in results.items())
            + f"; load {os.getloadavg()[0]:.1f}", flush=True)


if __name__ == "__main__":
    main()
