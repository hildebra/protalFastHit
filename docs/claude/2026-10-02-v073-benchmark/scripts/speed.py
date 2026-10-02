#!/usr/bin/env python3
"""speed.py [BENCH_DIR] [RESULTS ...] - wall time, CPU time and the load average before the run, per read type, depth
(scenario's pairs or bases) and version, from score.py's runs.tsv (default: BENCH_DIR/results_v073): the means over
both databases. Prints a markdown table per results folder.
"""
import collections
import csv
import os
import re
import statistics
import sys

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
for name in sys.argv[2:] or ["results_v073"]:
    rows = collections.defaultdict(list)
    with open(os.path.join(B, name, "runs.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            depth = re.search(r"_([pb]\d+)", r["scenario"]).group(1)
            rows[(r["reads"], depth, r["variant"])].append(r)
    variants = list(dict.fromkeys(k[2] for k in rows))
    print(f"## {name}: wall s / CPU s (load before), means over both databases\n")
    print("| reads | depth | runs | " + " | ".join(variants) + " |")
    print("|---|---|---|" + "---|" * len(variants))
    for reads, depth in sorted({k[:2] for k in rows}, key=lambda k: (k[0], int(k[1][1:]))):
        cells, n = [], 0
        for v in variants:
            rs = rows.get((reads, depth, v))
            if not rs:
                cells.append("-")
                continue
            n = len(rs)
            mean = lambda key: statistics.mean(float(r[key]) for r in rs if r[key] not in ("", "nan"))
            cells.append(f"{mean('wall_s'):.1f} / {mean('cpu_s'):.0f} ({mean('load_before'):.1f})")
        print(f"| {reads} | {depth} | {n} | " + " | ".join(cells) + " |")
    print()
