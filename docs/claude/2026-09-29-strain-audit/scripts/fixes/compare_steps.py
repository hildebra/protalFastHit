#!/usr/bin/env python3
"""Compare trees of the phylogeny audit's runs between protal/qcmsa versions.

usage: compare_steps.py TAG [TAG...]   (e.g. L M: reads results_L.jsonl and results_Lk04.jsonl, ...)
The baseline (7c6b2f8) is results.jsonl / results_k04.jsonl. For each tag: the raw MSA and the
qcmsa-filtered one, where the tag's runs have them.
"""
import json
import os
import sys

P = os.path.expanduser("~/audit5/phylo/runs")
RUNS = ["base20", "d2", "d2r2", "d2r3", "d3", "d3r2", "d5", "d10", "d50", "uneven", "mixed", "congener", "congener5"]


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
                continue  # Cferv's trees come from the --knob 0.4 runs
            out[(r["species"], r["variant"])] = r
    return out


def cell(r):
    if r is None:
        return "-"
    if r.get("status") != "ok":
        return f"{r.get('n_seqs', '?')} {r.get('status')}"
    return (f"{r['n_seqs']:>2} RF{r['rf']:<2} term{r['term_ratio']:.2f} pat{r['pat_ratio']:.2f} "
            f"T95 {r.get('n_true_splits_ge95')}/{r.get('n_true_splits')}")


def main():
    tags = sys.argv[1:]
    head = ["run", "sp", "raw 7c6b2f8", "filt 7c6b2f8"]
    for t in tags:
        head += [f"raw {t}", f"filt {t}"]
    print("\t".join(head))
    for run in RUNS:
        base = load(run, "")
        new = [load(run, t) for t in tags]
        for sp in ("Malpha", "Tone", "Cferv"):
            if (sp, "raw") not in base and not any((sp, "filt") in n for n in new):
                continue
            row = [run, sp, cell(base.get((sp, "raw"))), cell(base.get((sp, "filt")))]
            for n in new:
                row += [cell(n.get((sp, "raw"))), cell(n.get((sp, "filt")))]
            print("\t".join(row))


if __name__ == "__main__":
    main()
