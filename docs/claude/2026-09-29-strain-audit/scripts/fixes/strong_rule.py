#!/usr/bin/env python3
"""Counts of one strong-evidence rule on labelled dumps: flagged present-but-unreported (FN) and
flagged absent (TN, FP) taxa. usage: strong_rule.py DIR..."""
import glob
import os
import sys

RULES = {
    "d1 hgf.9 top.99": lambda r: r["depth"] >= 1 and r["hit_gene_fraction"] >= 0.9 and r["top_identity"] >= 0.99,
    "d1 hgf.9 top.98": lambda r: r["depth"] >= 1 and r["hit_gene_fraction"] >= 0.9 and r["top_identity"] >= 0.98,
    "d1 hgf.9 top.98 lis.5": lambda r: r["depth"] >= 1 and r["hit_gene_fraction"] >= 0.9 and r["top_identity"] >= 0.98
                                       and r["low_identity_share"] <= 0.5,
    "d1 hgf.8 top.98 lis.5": lambda r: r["depth"] >= 1 and r["hit_gene_fraction"] >= 0.8 and r["top_identity"] >= 0.98
                                       and r["low_identity_share"] <= 0.5,
    "d1 hgf.5 top.98 lis.5": lambda r: r["depth"] >= 1 and r["hit_gene_fraction"] >= 0.5 and r["top_identity"] >= 0.98
                                       and r["low_identity_share"] <= 0.5,
}
KEYS = ["depth", "hit_gene_fraction", "top_identity", "low_identity_share", "probability"]
for d in sys.argv[1:]:
    rows = []
    for path in glob.glob(os.path.join(d, "*.truth_annotated")):
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                rows.append({"truth": r["truth"], "prediction": r["prediction"], **{k: float(r[k]) for k in KEYS}})
    below = [r for r in rows if r["prediction"] == "0"]
    fn = [r for r in below if r["truth"] == "1"]
    tn = [r for r in below if r["truth"] == "0"]
    print(f"== {d}: {len(fn)} present below the knob, {len(tn)} absent below it")
    for name, rule in RULES.items():
        print(f"  {name:24s} FN {sum(map(rule, fn)):3d}  TN {sum(map(rule, tn)):2d}")
