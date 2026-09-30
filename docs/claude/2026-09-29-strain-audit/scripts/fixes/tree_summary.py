#!/usr/bin/env python3
"""Per step, per depth group and MSA variant: wrong trees, true splits with >= 95 bootstrap,
terminal and patristic length ratios, and sequences qcmsa removed.

usage: tree_summary.py TAG...   ("" or base = the audit's baseline 7c6b2f8, results.jsonl)
"""
import json
import os
import sys

P = os.path.expanduser("~/audit5/phylo/runs")
GROUPS = [("2x", ["d2", "d2r2", "d2r3"]), ("3x", ["d3", "d3r2"]),
          ("5-50x", ["d5", "d10", "d50", "base20", "uneven"]), ("mixed, congeners", ["mixed", "congener", "congener5"])]


def load(run, tag):
    out = {}
    for suffix in ("", "k04"):
        t = tag + suffix
        path = os.path.join(P, run, f"results_{t}.jsonl" if t else "results.jsonl")
        if not os.path.exists(path):
            continue
        for line in open(path):
            r = json.loads(line)
            if (r["species"] == "Cferv") != (suffix == "k04"):
                continue
            out[(r["species"], r["variant"])] = r
    return out


print("step\tdepth\tvariant\ttrees\twrong\tT95\tterminal\tpatristic\tseqs removed")
for tag in sys.argv[1:]:
    t = "" if tag == "base" else tag
    for name, runs in GROUPS:
        for variant in ("raw", "filt"):
            n = wrong = t95 = splits = removed = 0
            term, pat = [], []
            for run in runs:
                recs = load(run, t)
                for (sp, v), r in recs.items():
                    if v != variant or r.get("status") != "ok":
                        continue
                    n += 1
                    wrong += r["rf"] > 0
                    t95 += r.get("n_true_splits_ge95") or 0
                    splits += r.get("n_true_splits") or 0
                    term.append(r["term_ratio"])
                    pat.append(r["pat_ratio"])
                    if variant == "filt" and (sp, "raw") in recs:
                        removed += recs[(sp, "raw")]["n_seqs"] - r["n_seqs"]
            if n == 0:
                continue
            print(f"{tag}\t{name}\t{variant}\t{n}\t{wrong}\t{t95}/{splits}\t{min(term):.2f}-{max(term):.2f}\t"
                  f"{min(pat):.2f}-{max(pat):.2f}\t{removed if variant == 'filt' else ''}")
