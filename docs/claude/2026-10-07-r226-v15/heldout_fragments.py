#!/usr/bin/env python3
"""How many of the fragments error_reads.py tracked in a sample came from held-out species? error_reads.py keeps a
Python object for every fragment that has a reason (on an FP/FN taxon, seeded on an FN/unseen taxon, or drawn from an
FN/unseen species' genome); in v15 the held-out species counted as unseen, so all their reads with a record were
tracked. Per sample: the summary's fragments, and the own fragments of the unseen species that are held out (an
upper bound of what --heldout removes: a fragment can have another reason too).

    python3 heldout_fragments.py ~/v15/model_logs/error_reads/ont --heldout local/v15/heldout_species.txt
"""
import argparse
import os

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("dir")
    ap.add_argument("--heldout", required=True)
    opts = ap.parse_args()
    with open(opts.heldout) as fh:
        heldout = {line.split("\t")[0].strip() for line in fh if line.strip()}
    summary = pd.read_csv(os.path.join(opts.dir, "summary.tsv"), sep="\t")
    rows = []
    for _, s in summary.iterrows():
        path = os.path.join(opts.dir, s["set"], s["point"], s["sample"] + ".taxa.tsv")
        if not os.path.exists(path):
            continue
        t = pd.read_csv(path, sep="\t", usecols=["error", "taxon_name", "own_fragments"])
        t["own_fragments"] = pd.to_numeric(t["own_fragments"], errors="coerce").fillna(0)
        held = t[(t["error"] == "unseen") & t["taxon_name"].isin(heldout)]["own_fragments"].sum()
        rows.append({"set": s["set"], "scenario": s["scenario"], "sample": s["sample"], "fragments": s["fragments"],
                     "held-out species' own": int(held), "source MB": s["source_sam_bytes"] / 1e6})
    d = pd.DataFrame(rows)
    g = d.groupby(["set", "scenario"]).agg(samples=("sample", "size"), fragments=("fragments", "mean"),
                                           heldout=("held-out species' own", "mean"), source_mb=("source MB", "mean"))
    g["held-out share"] = (g["heldout"] / g["fragments"]).round(3)
    g["fragments left"] = (g["fragments"] - g["heldout"]).round(0)
    pd.set_option("display.width", 200)
    print(f"# {opts.dir}: per sample, means")
    print(g.round(1).to_string())


if __name__ == "__main__":
    main()
