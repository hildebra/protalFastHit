#!/usr/bin/env python3
"""Trees of several MSA variants of one tag, per run and species.
usage: variants_cmp.py TAG VARIANT [VARIANT...]   (reads runs/<run>/results_TAG[k04].jsonl)"""
import json
import os
import sys

P = os.path.expanduser("~/audit5/phylo/runs")
RUNS = ["base20", "d2", "d2r2", "d2r3", "d3", "d3r2", "d5", "d10", "d50", "uneven", "mixed", "congener", "congener5"]
tag, variants = sys.argv[1], sys.argv[2:]
print("\t".join(["run", "sp"] + variants))
for run in RUNS:
    recs = {}
    for suffix in ("", "k04"):
        path = os.path.join(P, run, f"results_{tag}{suffix}.jsonl")
        if not os.path.exists(path):
            continue
        for line in open(path):
            r = json.loads(line)
            if (r["species"] == "Cferv") != (suffix == "k04"):
                continue
            recs[(r["species"], r["variant"])] = r
    for sp in ("Malpha", "Tone", "Cferv"):
        cells = []
        for v in variants:
            r = recs.get((sp, v))
            if r is None:
                cells.append("-")
            elif r.get("status") != "ok":
                cells.append(f"{r.get('n_seqs', '?')} {r.get('status')}")
            else:
                cells.append(f"{r['n_seqs']:>2} RF{r['rf']:<2} term{r['term_ratio']:.2f} pat{r['pat_ratio']:.2f} "
                             f"T95 {r.get('n_true_splits_ge95')}/{r.get('n_true_splits')}")
        if any(c != "-" for c in cells):
            print("\t".join([run, sp] + cells))
