#!/usr/bin/env python3
"""From evidence.py's calls.tsv: the calls of each run, set and read type by their fragments (TP and FP per bin), and the
median F1 per run, set and read type if the calls with fewer than N fragments were dropped (abundance and closed
truths; a dropped TP becomes an FN), or with a share of the profile below S.

    python3 bins.py OUT_DIR/calls.tsv OUT_DIR        writes bins.tsv, min_fragments.tsv and min_share.tsv
"""
import collections
import csv
import os
import statistics
import sys

BINS = [(1, 2), (3, 5), (6, 10), (11, 25), (26, 100), (101, 1000), (1001, 10 ** 12)]
MINIMA = [0, 3, 6, 11, 26, 51, 101]


def main(calls, out):
    rows = [r for r in csv.DictReader(open(calls, newline=""), delimiter="\t") if r["run"] != "0.8.5_956e18b_v0"]
    bins = collections.Counter()
    samples = collections.defaultdict(list)
    for r in rows:
        samples[(r["run"], r["set"], r["read_type"], r["sample"])].append(r)
        if r["status"] in ("TP", "FP") and r["fragments"] != "NA":
            n = int(r["fragments"])
            label = next(f"{a}-{b}" if b < 10 ** 12 else f">{a - 1}" for a, b in BINS if a <= n <= b)
            bins[(r["run"], r["set"], r["read_type"], r["status"], label)] += 1
    with open(os.path.join(out, "bins.tsv"), "w") as fh:
        labels = [f"{a}-{b}" if b < 10 ** 12 else f">{a - 1}" for a, b in BINS]
        fh.write("\t".join(["run", "set", "read_type", "status"] + labels) + "\n")
        for key in sorted({k[:4] for k in bins}):
            fh.write("\t".join(list(key) + [str(bins[key + (l,)]) for l in labels]) + "\n")
    f1 = collections.defaultdict(list)
    for (run, set_name, read_type, sample), members in samples.items():
        statuses = {m["status"] for m in members}
        if not statuses <= {"TP", "FP", "FN"} or "NA" in {m["fragments"] for m in members if m["status"] != "FN"}:
            continue  # dropout/negative samples, or calls without evidence (sylph, 0.6.0a's cut logs)
        for minimum in MINIMA:
            tp = sum(m["status"] == "TP" and int(m["fragments"]) >= minimum for m in members)
            fp = sum(m["status"] == "FP" and int(m["fragments"]) >= minimum for m in members)
            fn = sum(m["status"] == "FN" for m in members) + sum(m["status"] == "TP" and int(m["fragments"]) < minimum
                                                                for m in members)
            f1[(run, set_name, read_type, minimum)].append(2 * tp / (2 * tp + fp + fn) if tp + fp + fn else 1.0)
    with open(os.path.join(out, "min_fragments.tsv"), "w") as fh:
        fh.write("\t".join(["run", "set", "read_type", "samples"] + [f"f1_min{m}" for m in MINIMA]) + "\n")
        for key in sorted({k[:3] for k in f1}):
            fh.write("\t".join(list(key) + [str(len(f1[key + (0,)]))] +
                               [f"{statistics.median(f1[key + (m,)]):.4f}" for m in MINIMA]) + "\n")


SHARES = [0, 1e-6, 3e-6, 1e-5, 3e-5, 1e-4]


def share_floor(calls, out):
    """The median F1 per run, set and read type if the calls below a share of the profile (profile_share: the placed
    share, as scored) were dropped; written to min_share.tsv. Every run, sylph's too (no evidence needed)."""
    samples = collections.defaultdict(list)
    for r in csv.DictReader(open(calls, newline=""), delimiter="\t"):
        if r["run"] != "0.8.5_956e18b_v0":
            samples[(r["run"], r["set"], r["read_type"], r["sample"])].append(r)
    f1 = collections.defaultdict(list)
    for (run, set_name, read_type, sample), members in samples.items():
        if not {m["status"] for m in members} <= {"TP", "FP", "FN"}:
            continue
        for floor in SHARES:
            kept = [m for m in members if m["status"] != "FN" and float(m["profile_share"]) >= floor]
            tp = sum(m["status"] == "TP" for m in kept)
            fp = sum(m["status"] == "FP" for m in kept)
            fn = sum(m["status"] != "FP" for m in members) - tp
            f1[(run, set_name, read_type, floor)].append(2 * tp / (2 * tp + fp + fn) if tp + fp + fn else 1.0)
    with open(os.path.join(out, "min_share.tsv"), "w") as fh:
        fh.write("\t".join(["run", "set", "read_type", "samples"] + [f"f1_share{s:g}" for s in SHARES]) + "\n")
        for key in sorted({k[:3] for k in f1}):
            fh.write("\t".join(list(key) + [str(len(f1[key + (0,)]))] +
                               [f"{statistics.median(f1[key + (s,)]):.4f}" for s in SHARES]) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:3])
    share_floor(*sys.argv[1:3])
