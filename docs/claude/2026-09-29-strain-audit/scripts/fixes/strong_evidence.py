#!/usr/bin/env python3
"""Which rule flags species that are present but score below the knob, and nothing absent?

usage: strong_evidence.py DIR [DIR...]   (directories of <sample>.profile.truth_annotated)
For the taxa below the knob (prediction 0): how many of the present (truth 1) and absent (truth 0)
ones each candidate rule flags, and feature quantiles of both groups.
"""
import glob
import os
import sys

rows = []
for d in sys.argv[1:]:
    for path in glob.glob(os.path.join(d, "*.truth_annotated")):
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                r["_sample"] = os.path.basename(path)
                rows.append(r)

num = lambda r, k: float(r[k])
fn = [r for r in rows if r["truth"] == "1" and r["prediction"] == "0"]
tn = [r for r in rows if r["truth"] == "0" and r["prediction"] == "0"]
tp = [r for r in rows if r["truth"] == "1" and r["prediction"] == "1"]
fp = [r for r in rows if r["truth"] == "0" and r["prediction"] == "1"]
print(f"rows {len(rows)}: TP {len(tp)} FN {len(fn)} FP {len(fp)} TN {len(tn)}")

feats = ["depth", "hit_gene_fraction", "gene_presence_ratio", "top_identity", "identity", "low_identity_share",
         "probability", "total_hits"]


def q(vals, p):
    s = sorted(vals)
    return s[min(len(s) - 1, int(p * (len(s) - 1) + 0.5))] if s else float("nan")


for name, group in (("FN (present, below knob)", fn), ("TN (absent, below knob)", tn)):
    print(f"\n{name}: {len(group)}")
    for f in feats:
        v = [num(r, f) for r in group]
        print(f"  {f:22s} " + " ".join(f"{q(v, p):9.4g}" for p in (0.05, 0.25, 0.5, 0.75, 0.95)))

rules = []
for depth in (0.5, 1, 2):
    for hgf in (0.3, 0.5, 0.7, 0.9):
        for top in (0.95, 0.97, 0.98, 0.99):
            rules.append((depth, hgf, top))
print("\nrule (depth, hit_gene_fraction, top_identity): flagged FN / flagged TN")
best = []
for depth, hgf, top in rules:
    ok = lambda r: num(r, "depth") >= depth and num(r, "hit_gene_fraction") >= hgf and num(r, "top_identity") >= top
    a = sum(ok(r) for r in fn)
    b = sum(ok(r) for r in tn)
    best.append((b, -a, depth, hgf, top))
for b, a, depth, hgf, top in sorted(best)[:25]:
    print(f"  depth>={depth} hgf>={hgf} top>={top}: {-a}/{len(fn)} FN, {b}/{len(tn)} TN")
print("\nflagged TN under depth>=1 hgf>=0.5 top>=0.98:")
for r in tn:
    if num(r, "depth") >= 1 and num(r, "hit_gene_fraction") >= 0.5 and num(r, "top_identity") >= 0.98:
        print("  ", r["_sample"], r["taxon_name"], {f: r[f] for f in feats})
