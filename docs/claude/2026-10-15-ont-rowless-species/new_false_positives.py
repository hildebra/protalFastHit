#!/usr/bin/env python3
"""The false positives a lower floor adds (docs/claude/2026-10-15-ont-rowless-species): per set, the (sample, taxon)
calls of absent taxa in a build that the first build does not make, whether the first build had a row for them, and
whether the sample holds a congener left out of the training database (meta_novel_congener).

Usage: new_false_positives.py LABEL=BUILD_DIR ... (the first is the reference)
"""
import os
import sys

import pandas as pd

COLUMNS = ["meta_design", "meta_sample", "meta_novel_congener", "taxon", "truth", "set", "call"]


def main():
    builds = []
    for arg in sys.argv[1:]:
        label, build = arg.split("=", 1)
        c = pd.read_csv(os.path.join(build, "model_logs", "trained_model_ont.calls.tsv.gz"), sep="\t", usecols=COLUMNS)
        builds.append((label, c))
    (ref_label, ref), others = builds[0], builds[1:]
    for label, c in others:
        for which in ("training", "test"):
            r, d = ref[ref["set"] == which], c[c["set"] == which]
            key = lambda t: set(zip(t["meta_sample"], t["taxon"]))
            fp_ref, fp = key(r[(r["call"] == 1) & (r["truth"] == 0)]), d[(d["call"] == 1) & (d["truth"] == 0)]
            new = fp[[k not in fp_ref for k in zip(fp["meta_sample"], fp["taxon"])]]
            had_row = [k in key(r) for k in zip(new["meta_sample"], new["taxon"])]
            print(f"{which}, {label} against {ref_label}: false positives {len(fp_ref)} -> {len(fp)}; "
                  f"{len(new)} new ({sum(had_row)} with a row at {ref_label}), {len(fp_ref - key(fp))} gone; "
                  f"of the new, {int((new['meta_novel_congener'] == 1).sum())} in samples with a congener left out")


if __name__ == "__main__":
    main()
