#!/usr/bin/env python3
"""The share of a long-read sample's time that holds the GIL: a pure-Python counter on a second thread counts
alone for a few seconds, then beside one sample (long_read_templates + hifi_reads.simulate, or the templates
only); the counter's rate lost is about the share the sample held the GIL. With N threads in one process a
sample then takes at least N x that share of its time, so it bounds the thread pool's speed-up at 1 / share.

usage: gil_share.py MANIFEST OUT [--bases 30000000]
"""
import argparse
import csv
import os
import sys
import threading
import time

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
import collect_training_data as ctd  # noqa: E402
import hifi_reads  # noqa: E402


def count_while(stop):
    n = 0
    while not stop.is_set():
        for _ in range(1000):
            n += 1
    return n


def counter_rate(work):
    """Counts per second of a counter thread while work() runs (or for 5 s if work is None), and work's time."""
    stop, result = threading.Event(), {}
    t = threading.Thread(target=lambda: result.update(n=count_while(stop)))
    began = time.time()
    t.start()
    if work is None:
        time.sleep(5)
    else:
        work()
    seconds = time.time() - began
    stop.set()
    t.join()
    return result["n"] / seconds, seconds


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest")
    p.add_argument("out")
    p.add_argument("--bases", type=int, default=30_000_000)
    opts = p.parse_args()
    os.makedirs(opts.out, exist_ok=True)
    with open(opts.manifest) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    first = rows[0]["sample"]
    genomes = [{"genome": r["genome"], "fasta": r["fasta_path"],
                "weight": float(r["relative_abundance"]) * float(r["genome_length"])} for r in rows if r["sample"] == first]
    setup = ctd.parse_long_setup("hifi:15000:3000:3")
    task = {"seed": 5, "genomes": genomes, "setup": setup, "bases": opts.bases,
            "templates": os.path.join(opts.out, "t.fa"), "out": os.path.join(opts.out, "r.fq.gz")}
    alone, _ = counter_rate(None)
    print("work\twall_s\tGIL share\tthread pool speed-up at most", flush=True)
    for name, work in (("templates", lambda: ctd.long_read_templates(task)),
                       ("hifi_reads.simulate", lambda: hifi_reads.simulate(task["templates"], task["out"], 3.0, 5))):
        rate, seconds = counter_rate(work)
        share = max(0.0, 1 - rate / alone)
        print(f"{name}\t{seconds:.1f}\t{share:.2f}\t{(1 / share if share else float('inf')):.1f}", flush=True)


if __name__ == "__main__":
    main()
