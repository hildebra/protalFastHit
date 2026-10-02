#!/usr/bin/env python3
"""fragment_depth_summary.py OUT - fragment_depth_compare.sh's two arms side by side: per read setup, the calls (TP,
FP, FN against the truth), how many taxa changed their call, and the present taxa's depth new / old (median, and the
5th and 95th percentile), with the relative abundance of the present taxa (depth over the sample's passing depth)
new against old."""
import collections
import csv
import glob
import os
import statistics
import sys

csv.field_size_limit(1 << 30)
out = sys.argv[1]


def load(arm):
    taxa = {}
    for path in glob.glob(os.path.join(out, arm, "*.profile.truth_annotated")):
        sample = os.path.basename(path).split(".profile")[0]
        for r in csv.DictReader(open(path), delimiter="\t"):
            taxa[(sample, r["taxon"])] = (int(r["truth"]), int(r["prediction"]), float(r["depth"]))
    return taxa


old, new = load("old"), load("new")


def setup(sample):
    return sample.split("_p")[0]


def quantile(v, q):
    v = sorted(v)
    return v[min(len(v) - 1, int(q * len(v)))]


print("| reads (fragments) | samples | TP old / new | FP old / new | FN old / new | calls changed | present taxa's depth new / old: median (5-95%) | their abundance new / old: median (5-95%) |")
print("|---|---|---|---|---|---|---|---|")
frag = {"rl100": "2x100 (300 +- 40)", "rl150": "2x150 (350 +- 50)", "rl250": "2x250 (550 +- 50)"}
for s in sorted({setup(k[0]) for k in old}):
    keys = [k for k in old if setup(k[0]) == s and k in new]
    count = lambda taxa, t, p: sum(1 for k in keys if taxa[k][0] == t and taxa[k][1] == p)
    changed = sum(1 for k in keys if old[k][1] != new[k][1])
    ratios = [new[k][2] / old[k][2] for k in keys if old[k][0] == 1 and old[k][2] > 0]
    totals = collections.defaultdict(lambda: [0.0, 0.0])
    for k in keys:
        if old[k][0] == 1:
            totals[k[0]][0] += old[k][2]
            totals[k[0]][1] += new[k][2]
    abundance = [(new[k][2] / totals[k[0]][1]) / (old[k][2] / totals[k[0]][0]) for k in keys
                 if old[k][0] == 1 and old[k][2] > 0 and totals[k[0]][0] > 0]
    samples = len({k[0] for k in keys})
    fmt = lambda v: f"{statistics.median(v):.4f} ({quantile(v, 0.05):.4f}-{quantile(v, 0.95):.4f})" if v else "-"
    print(f"| {frag.get(s, s)} | {samples} | {count(old, 1, 1)} / {count(new, 1, 1)} | {count(old, 0, 1)} / {count(new, 0, 1)} | "
          f"{count(old, 1, 0)} / {count(new, 1, 0)} | {changed} | {fmt(ratios)} | {fmt(abundance)} |")
