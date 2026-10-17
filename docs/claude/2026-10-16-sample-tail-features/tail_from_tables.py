#!/usr/bin/env python3
"""The sample tail features (protal's context::SampleComplexityOf since 2026-10-16) computed from a collection's table:
sample_thin_share, the sample's taxa of 1-2 fragments over its taxa with fragments, and sample_thin_close_share, the
share of those thin taxa whose identity is within --margin (0.01, protal's context::kThinIdentityMargin) of
sample_identity. The table's `fragments` and `identity` are
the values protal summarises (Fragments(), BaseIdentity()), and its rows the sample's taxa: 10^sample_log_taxa is checked
against the rows with fragments.

Usage:
  tail_from_tables.py [--margin M] summary TABLE [TABLE ...]       per scenario sample: depth, taxa, the two features, its thin
                                                      present species (meta_marker_reads 1-2 among truth 1)
  tail_from_tables.py [--margin M] add TABLE OUT                   TABLE with the two columns added (after sample_identity)
"""
import argparse
import math
import sys

import numpy as np
import pandas as pd

THIN_FRAGMENTS = 2
THIN_IDENTITY_MARGIN = 0.01  # protal's context::kThinIdentityMargin
SIGMAS = (1.0, 1.5, 2.0, 2.5)


def tail_of(table, margin=THIN_IDENTITY_MARGIN):
    """{sample: (thin_share, thin_close_share, taxa with fragments, 10^sample_log_taxa)} of a table's rows."""
    out = {}
    for sample, g in table.groupby("meta_sample", sort=False):
        f = g["fragments"].to_numpy(dtype=float)
        ident = g["identity"].to_numpy(dtype=float)
        with_fragments = f >= 1
        thin = with_fragments & (f <= THIN_FRAGMENTS)
        sample_identity = float(g["sample_identity"].iloc[0])
        close = thin & (ident >= sample_identity - margin)
        n = int(with_fragments.sum())
        out[sample] = (thin.sum() / max(n, 1), close.sum() / thin.sum() if thin.sum() else 0.0, n,
                       10 ** float(g["sample_log_taxa"].iloc[0]))
    return out


def summary(paths, margin=THIN_IDENTITY_MARGIN):
    for path in paths:
        t = pd.read_csv(path, sep="\t", usecols=lambda c: c in {
            "meta_sample", "meta_scenario", "meta_marker_reads", "truth", "fragments", "identity", "sample_identity",
            "sample_log_taxa", "sample_log_fragments"}, dtype={"meta_sample": str, "meta_scenario": str},
            keep_default_na=False, na_values=[""])
        tails = tail_of(t, margin)
        bad = [s for s, v in tails.items() if abs(v[2] - v[3]) > 0.5]
        print(f"## {path}: {len(tails)} samples; taxa with fragments == 10^sample_log_taxa in "
              f"{len(tails) - len(bad)} of them" + (f" (not: {', '.join(bad[:5])})" if bad else ""))
        print("\t".join(["sample", "sigma", "sample_log_fragments", "sample_log_taxa", "sample_thin_share",
                         "sample_thin_close_share", "thin present (1-2 marker reads)", "absent rows"]))
        scen = t[t["meta_scenario"].astype(str).isin(["soil", "soil_shallow", "moderate"])]
        for sample, g in sorted(scen.groupby("meta_sample"), key=lambda x: (x[0].rsplit("_s_", 1)[0],
                                                                         int(x[0].rsplit("_s_", 1)[1]))):
            n = int(sample.rsplit("_s_", 1)[1])
            mr = pd.to_numeric(g["meta_marker_reads"], errors="coerce")
            thin_present = int(((g["truth"] == 1) & (mr <= 2)).sum())
            v = tails[sample]
            print("\t".join([sample, f"{SIGMAS[(n - 1) % 4]:g}", f"{float(g['sample_log_fragments'].iloc[0]):.3f}",
                             f"{float(g['sample_log_taxa'].iloc[0]):.3f}", f"{v[0]:.3f}", f"{v[1]:.3f}", str(thin_present),
                             str(int((g['truth'] == 0).sum()))]))
        print()


def add(path, out, margin=THIN_IDENTITY_MARGIN):
    t = pd.read_csv(path, sep="\t", dtype={"meta_sample": str}, keep_default_na=False, na_values=[""], low_memory=False)
    tails = tail_of(t, margin)
    thin = t["meta_sample"].map(lambda s: tails[s][0])
    close = t["meta_sample"].map(lambda s: tails[s][1])
    at = list(t.columns).index("sample_identity") + 1
    t.insert(at, "sample_thin_close_share", close)
    t.insert(at, "sample_thin_share", thin)
    t.to_csv(out, sep="\t", index=False, float_format="%.9g")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="command", required=True)
    s = sub.add_parser("summary")
    s.add_argument("tables", nargs="+")
    a = sub.add_parser("add")
    a.add_argument("table")
    a.add_argument("out")
    ap.add_argument("--margin", type=float, default=THIN_IDENTITY_MARGIN,
                    help=f"the identity margin of a close thin taxon (default {THIN_IDENTITY_MARGIN}, protal's)")
    opts = ap.parse_args()
    if opts.command == "summary":
        summary(opts.tables, opts.margin)
    else:
        add(opts.table, opts.out, opts.margin)


if __name__ == "__main__":
    main()
