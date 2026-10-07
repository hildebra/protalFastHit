#!/usr/bin/env python3
"""Where the reads of the missed species went: error_reads.py's <sample>.taxa.tsv tables, summed by set, scenario and
error (FN: present, scored, not called; unseen: present in the sample and in the training database, no row at all).

A missed species is put in the first class it fits, by its own fragments (reads drawn from its genomes that left a
record in the SAM, error_reads.py's own_* columns) and the fragments that seeded on it:
  on it, MAPQ >= 4       a fragment's best record is on the species at MAPQ 4 or more (as the profiler counts them)
  on it, MAPQ < 4        best records on the species only below MAPQ 4 (shared with relatives)
  none                   no fragment of its own left a record (no read seeded on any marker gene)
  elsewhere              some of its fragments' best records are on other taxa (own_best_elsewhere)
  unaligned, seeded on it
                         its fragments aligned nowhere (own_unaligned), and fragments seeded on its genes and failed
                         to align to them (seeded_not_aligned, anyone's fragments by the ZF tag: most likely its own)
  unaligned elsewhere    its fragments aligned nowhere and none seeded on its genes: a long read from outside the
                         marker genes leaves an unmapped record from a few chance k-mers elsewhere
FP rows are counted, with the fragments on them, for scale.

error_reads.py takes "the training database's species" from training_db/genome2tiid.tsv, which the converter copies
unchanged (--from_db keeps every species in the taxonomy files and leaves the held-out ones' genes out), so its unseen
species include the held-out ones: novel to the training database, their reads can only land on relatives or nowhere.
With --heldout (the build's heldout_species.txt) those are counted apart as "unseen (held out)".

    python3 error_taxa_summary.py ~/v15/model_logs/error_reads/ont --heldout local/v15/heldout_species.txt
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

NUM = ["own_fragments", "own_best_on_taxon", "own_best_on_taxon_mapq4", "own_best_elsewhere", "own_unaligned",
       "fragments_on_taxon", "seeded_not_aligned"]


def scenario_of(point):
    for sc in ("soil_shallow", "soil", "gut", "host"):
        if point.startswith(f"sc_{sc}_"):
            return sc
    return "design"


def classify(t):
    own = t["own_fragments"]
    return np.select([t["own_best_on_taxon_mapq4"] > 0, t["own_best_on_taxon"] > 0, own <= 0,
                      t["own_best_elsewhere"] > 0, t["seeded_not_aligned"] > 0],
                     ["on it, MAPQ >= 4", "on it, MAPQ < 4", "none", "elsewhere", "unaligned, seeded on it"],
                     "unaligned elsewhere")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("dir", help="model_logs/error_reads/<read type>")
    ap.add_argument("--heldout", help="the build's heldout_species.txt: unseen held-out species counted apart")
    opts = ap.parse_args()
    heldout = set()
    if opts.heldout:
        with open(opts.heldout) as fh:
            heldout = {line.split("\t")[0].strip() for line in fh if line.strip()}
    frames = []
    for path in glob.glob(os.path.join(opts.dir, "*", "*", "*.taxa.tsv")):
        t = pd.read_csv(path, sep="\t", dtype={"taxid": str})
        if t.empty:
            continue
        t["point"] = os.path.basename(os.path.dirname(path))
        frames.append(t)
    t = pd.concat(frames, ignore_index=True)
    for c in NUM:
        t[c] = pd.to_numeric(t[c], errors="coerce").fillna(0)
    t["scenario"] = t["point"].map(scenario_of)
    t.loc[(t["error"] == "unseen") & t["taxon_name"].isin(heldout), "error"] = "unseen (held out)"
    samples = t.groupby(["set", "scenario"])["sample"].nunique().rename("samples")
    missed = t[t["error"].isin(["FN", "unseen", "unseen (held out)"])].copy()
    missed["class"] = classify(missed)
    tab = missed.groupby(["set", "scenario", "error", "class"]).size().unstack("class", fill_value=0)
    tab = tab.join(samples, on=["set", "scenario"])
    per = tab.drop(columns="samples").div(tab["samples"], axis=0).round(1)
    per.insert(0, "samples", tab["samples"])
    per["all"] = tab.drop(columns="samples").sum(axis=1).div(tab["samples"]).round(1)
    pd.set_option("display.width", 250)
    print(f"# {opts.dir}: {t['sample'].nunique()} samples; missed species per sample by where their own reads went")
    print(per.to_string())
    fp = t[t["error"] == "FP"].groupby(["set", "scenario"]).agg(FP=("taxid", "size"), fragments=("fragments_on_taxon", "median"))
    fp = fp.join(samples)
    fp["FP per sample"] = (fp["FP"] / fp["samples"]).round(1)
    print("\n# FP per sample (median fragments on an FP taxon)")
    print(fp.to_string())
    # the unseen species whose reads are on it at MAPQ >= 4 had reads a row could be made from
    u = missed[(missed["error"] == "unseen") & (missed["class"] == "on it, MAPQ >= 4")]
    if len(u):
        print("\n# unseen species with fragments best on them at MAPQ >= 4: those fragments (quantiles 50/90/99%)")
        print(u.groupby(["set", "scenario"])["own_best_on_taxon_mapq4"].quantile([0.5, 0.9, 0.99]).unstack().to_string())


if __name__ == "__main__":
    main()
