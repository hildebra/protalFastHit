#!/usr/bin/env python3
"""stride_compare.py - the hash sample of a read's cores against the stride pick, from a protal run with --stage1_probe.

Reads <dir>/<sample>_stage1_probe.tsv (the rules of the hash sample: the cores whose hash falls in 1/c of its range) and
<dir>/<sample>_stage1_probe_stride.tsv (the same rules for every c-th core of each mate in read order; protal since
2026-10-13), and with --profile the taxa table <dir>/<sample>_stage1_probe_taxa.tsv and the run's profile (the called
species). <dir> is a run's misc/ folder or measure_performance.sh's OUT_DIR/stage1_probe/. Writes Markdown:

  1. per c, for T = 14 (--score) and m = 1 and 2: reads kept, recall (the share of the reads with records kept) and
     lookups left (the share of the run's lookups a gate with the rule makes), hash against stride;
  2. the rule with the fewest lookups left at a recall of at least 99%, 99.5%, 99.8% and 99.9%, for each pick: the
     comparison at equal cost;
  3. with --profile: per c (T, m = 1), the called species keeping less than 90% of their reads with records.

  python3 stride_compare.py --dir OUT/stage1_probe --sample pe1 [--profile OUT/stage1_probe/pe1.plain.profile]

docs/claude/2026-10-13-stride-pick; the tallies: src/Core/Stage1Probe.h.
"""

import argparse
import csv
import os
import sys

LEVELS = [1, 2, 4, 8, 16, 32]
RECALLS = [0.99, 0.995, 0.998, 0.999]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def called_species(profile):
    """The lineages of the species a profile calls (its second column), without the '?' line."""
    called = set()
    with open(profile) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or line.startswith("#") or fields[0] == "?":
                continue
            called.add(fields[1])
    return called


def pct(value):
    return f"{100 * value:.2f}%"


def left(rule):
    """The share of lookups left: the table's (protal since 2026-10-13), else the estimate 1/c + kept x (1 - 1/c)."""
    if rule.get("lookups_left"):
        return float(rule["lookups_left"])
    c = int(rule["c"])
    return 1 / c + float(rule["share_kept"]) * (1 - 1 / c)


def keyed(rules):
    return {(int(r["c"]), int(r["min_score"]), int(r["min_hits"])): r for r in rules}


def cheapest(rules, recall):
    """The rule with the fewest lookups left at a recall of at least `recall` (ties: fewer reads kept)."""
    fitting = [r for r in rules if float(r["recall"]) >= recall]
    return min(fitting, key=lambda r: (left(r), float(r["share_kept"]))) if fitting else None


def losing_species(taxa, called, column, below=0.9):
    """The called species keeping less than `below` of their reads with records under the taxa table's `column`."""
    out = []
    for t in taxa:
        if t["name"] not in called or column not in t:
            continue
        aligned, kept = int(t["aligned"]), int(t[column])
        if aligned and kept < below * aligned:
            out.append((kept / aligned, aligned, t["name"]))
    return sorted(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dir", required=True, help="the run's misc folder, or measure_performance.sh's stage1_probe folder")
    parser.add_argument("--sample", required=True, help="the sample's name (the tables' prefix)")
    parser.add_argument("--profile", help="the sample's profile (for the called species)")
    parser.add_argument("--score", type=int, default=14, help="T of the per-c table and the species losses (14)")
    args = parser.parse_args()

    stride_path = os.path.join(args.dir, f"{args.sample}_stage1_probe_stride.tsv")
    if not os.path.exists(stride_path):
        print(f"{stride_path} is missing: the run's protal is older than the stride rules (2026-10-13)", file=sys.stderr)
        return 1
    hash_rules = read_tsv(os.path.join(args.dir, f"{args.sample}_stage1_probe.tsv"))
    stride_rules = read_tsv(stride_path)
    by_hash, by_stride = keyed(hash_rules), keyed(stride_rules)
    reads, aligned = int(hash_rules[0]["reads"]), int(hash_rules[0]["aligned"])

    print(f"# Hash sample against stride pick: {args.sample}\n")
    print(f"{reads:,} reads (pairs as one), {aligned:,} with records.\n")
    for m in (1, 2):
        print(f"## T = {args.score}, m = {m}: reads kept / recall / lookups left\n")
        print("| c | hash: kept | recall | lookups left | stride: kept | recall | lookups left | records lost, hash → stride |")
        print("|---|---|---|---|---|---|---|---|")
        for c in LEVELS:
            h, s = by_hash.get((c, args.score, m)), by_stride.get((c, args.score, m))
            if not h or not s:
                continue
            lost_h, lost_s = 1 - float(h["recall"]), 1 - float(s["recall"])
            print(f"| {c} | {pct(float(h['share_kept']))} | {pct(float(h['recall']))} | {pct(left(h))} | "
                  f"{pct(float(s['share_kept']))} | {pct(float(s['recall']))} | {pct(left(s))} | {pct(lost_h)} → {pct(lost_s)} |")
        print()

    print("## The fewest lookups at a recall\n")
    print("Over all rules (c, T, m) of each pick.\n")
    print("| recall at least | pick | c | T | m | reads kept | recall | lookups left |")
    print("|---|---|---|---|---|---|---|---|")
    for target in RECALLS:
        for name, rules in (("hash", hash_rules), ("stride", stride_rules)):
            r = cheapest(rules, target)
            if r is None:
                print(f"| {pct(target)} | {name} | none | | | | | |")
                continue
            print(f"| {pct(target)} | {name} | {r['c']} | {r['min_score']} | {r['min_hits']} | {pct(float(r['share_kept']))} | "
                  f"{pct(float(r['recall']))} | {pct(left(r))} |")
    print()

    taxa_path = os.path.join(args.dir, f"{args.sample}_stage1_probe_taxa.tsv")
    if args.profile and os.path.exists(taxa_path):
        called = called_species(args.profile)
        taxa = read_tsv(taxa_path)
        print(f"## Called species keeping less than 90% of their reads with records (T = {args.score}, m = 1)\n")
        print(f"{len(called)} species called.\n")
        print("| c | hash | stride |")
        print("|---|---|---|")
        for c in LEVELS:
            h = losing_species(taxa, called, f"c{c}_t{args.score}_m1")
            s = losing_species(taxa, called, f"s{c}_t{args.score}_m1")
            print(f"| {c} | {len(h)} | {len(s)} |")
        print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
