#!/usr/bin/env python3
"""How much two protal runs' profiles differ: calls, model probabilities and abundances.

usage: profile_diff.py DIR_A DIR_B   (each with profiles/<sample>.profile.log)
"""
import glob
import os
import sys


def load(d):
    out = {}
    for path in glob.glob(os.path.join(d, "profiles", "*.profile.log")):
        sample = os.path.basename(path)[:-len(".profile.log")]
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                out[(sample, r["Name"])] = (int(r["Predicted"]), float(r["Probability"]), float(r["Abundance"]))
    return out


a, b = load(sys.argv[1]), load(sys.argv[2])
keys = sorted(set(a) & set(b))
flips = [k for k in keys if a[k][0] != b[k][0]]
dp = [abs(a[k][1] - b[k][1]) for k in keys]
da = [abs(a[k][2] - b[k][2]) for k in keys if a[k][0] and b[k][0]]
ra = [abs(a[k][2] - b[k][2]) / a[k][2] for k in keys if a[k][0] and b[k][0] and a[k][2] > 0]
print(f"taxon rows: {len(keys)} (only in A: {len(set(a) - set(b))}, only in B: {len(set(b) - set(a))})")
print(f"presence calls that differ: {len(flips)}")
for k in flips[:20]:
    print(f"  {k[0]} {k[1]}: {a[k]} -> {b[k]}")
print(f"probability |diff|: mean {sum(dp) / len(dp):.4f}, max {max(dp):.4f}, rows > 0.01: {sum(d > 0.01 for d in dp)}")
if ra:
    print(f"abundance of taxa called in both: mean |rel diff| {sum(ra) / len(ra):.4f}, max {max(ra):.4f}, max |abs diff| {max(da):.5f}")
