#!/usr/bin/env python3
"""Where a long-read sample's Python time goes (cProfile, one PacBio HiFi sample of a sim_cost.sh community:
long_read_templates, then hifi_reads.simulate), by function: zlib and file reads, numpy's large array operations
release the GIL; bytes methods, list building and the interpreter itself hold it.

usage: gil_profile.py MANIFEST OUT [--bases 30000000] [--top 25]
"""
import argparse
import cProfile
import csv
import os
import pstats
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import collect_training_data as ctd  # noqa: E402
import hifi_reads  # noqa: E402

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("manifest")
p.add_argument("out")
p.add_argument("--bases", type=int, default=30_000_000)
p.add_argument("--top", type=int, default=25)
opts = p.parse_args()
os.makedirs(opts.out, exist_ok=True)
with open(opts.manifest) as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))
genomes = [{"genome": r["genome"], "fasta": r["fasta_path"],
            "weight": float(r["relative_abundance"]) * float(r["genome_length"])} for r in rows if r["sample"] == rows[0]["sample"]]
task = {"seed": 5, "genomes": genomes, "setup": ctd.parse_long_setup("hifi:15000:3000:3"), "bases": opts.bases,
        "templates": os.path.join(opts.out, "t.fa"), "out": os.path.join(opts.out, "r.fq.gz")}
for name, work in (("templates", lambda: ctd.long_read_templates(task)),
                   ("hifi_reads.simulate", lambda: hifi_reads.simulate(task["templates"], task["out"], 3.0, 5))):
    prof = cProfile.Profile()
    prof.runcall(work)
    stats = pstats.Stats(prof)
    print(f"== {name}: {stats.total_tt:.2f} s")
    stats.sort_stats("tottime").print_stats(opts.top)
