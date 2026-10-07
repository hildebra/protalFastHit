#!/usr/bin/env python3
"""Where a strain and a sister species should differ at the same marker identity: the mutation spectrum (the share of
mismatches at third codon positions, third_position_share) and the gene-rate pattern (excess_conserved_fast_ratio,
gene_divergence_dispersion), in the v15 tables by class of row (present representative / real strain / in-silico
strain; absent novel congener / other absent), at the knob's errors, by fragments; and whether the in-silico strains
show a real strain's spectrum. Also, from the fragment tables, whether the per-taxon divergence rescaled by the genes'
within- or between-species factors separates FP from FN taxa better than the raw divergence.

    python3 spectrum.py --records ~/v15/err_analysis --build local/v15 > spectrum.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import calls
import signatures

SPECTRUM = ["third_position_share", "excess_conserved_fast_ratio", "gene_divergence_dispersion",
            "conserved_fast_depth_ratio", "conserved_hit_share", "variant_sites_per_kb", "multiallelic_sites_per_kb",
            "failed_gene_share", "mate_lost_share", "identity", "excess_scaled_median", "low_mapq_share",
            "cluster_genomes_log10", "em_own_share", "lu_per_kb"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--types", default="pe,se,pb,ont")
    return p.parse_args(argv)


def by_class(t, cols, min_fragments):
    sub = t[t["fragments"] >= min_fragments]
    out = sub.groupby("class")[cols].median().T
    out["rows"] = ""
    n = sub["class"].value_counts()
    print(f"medians by class, rows with >= {min_fragments} fragments ({n.to_dict()}):")
    print(out.drop(columns="rows").round(4).to_string())


def aucs(t, cols, min_fragments):
    sub = t[t["fragments"] >= min_fragments]
    rows = []
    for col in cols:
        if col not in sub:
            continue
        r = {"feature": col}
        a = sub["error"].isin(["FN", "FP"])
        r["FN vs FP (FN high)"] = signatures.auc((sub.loc[a, "error"] == "FN").astype(int), sub.loc[a, col])
        b = sub["class"].isin(["present strain", "absent novel congener"])
        r["strain vs novel congener, all rows (strain high)"] = signatures.auc(
            (sub.loc[b, "class"] == "present strain").astype(int), sub.loc[b, col])
        band = b & sub["identity"].between(0.95, 0.985)
        r["same, identity 0.95-0.985"] = signatures.auc((sub.loc[band, "class"] == "present strain").astype(int),
                                                       sub.loc[band, col])
        r["n band"] = int(band.sum())
        rows.append(r)
    print(f"AUC, rows with >= {min_fragments} fragments:")
    print(pd.DataFrame(rows).round(3).to_string(index=False))


def insilico(t):
    sub = t[(t["class"].isin(["present strain", "present insilico"])) & (t["fragments"] >= 10)].copy()
    sub["band"] = pd.cut(sub["identity"], [0, 0.95, 0.97, 0.98, 0.99, 1.001])
    g = sub.groupby(["band", "class"], observed=True)["third_position_share"].agg(["size", "median"])
    print("third_position_share of present real and in-silico strains by identity band (>= 10 fragments):")
    print(g.round(3).to_string())


def rescaled(records, rt):
    f = pd.read_pickle(os.path.join(records, f"{rt}.fragments.pkl.gz"))
    c = f[f["counted"] & f["role"].isin(["FP", "FN"]) & f["gene"].notna()].copy()
    c = c[(c["role"] == "FP") | (c["relation"] == "own")]
    c["div"] = 1 - c["identity"]
    c["div_within"] = c["div"] / c["gene_within"].clip(lower=0.1)
    c["div_between"] = c["div"] / c["gene_between"].clip(lower=0.1)
    g = c.groupby(["sample", "taxon", "role"]).agg(n=("div", "size"), div=("div", "mean"), div_within=("div_within", "mean"),
                                                  div_between=("div_between", "mean")).reset_index()
    g["truth"] = (g["role"] == "FN").astype(int)
    for n in (1, 3, 10):
        s = g[g["n"] >= n]
        print(f"  taxa with >= {n} counted fragments ({len(s)}, FN {int(s['truth'].sum())}): AUC FN high: raw divergence "
              f"{signatures.auc(s['truth'], s['div']):.3f}, / within factor {signatures.auc(s['truth'], s['div_within']):.3f}, "
              f"/ between factor {signatures.auc(s['truth'], s['div_between']):.3f}")


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    for rt in opts.types.split(","):
        t, _ = calls.load(opts.build, rt)
        knob = calls.knob_of(opts.records, rt)
        t["class"] = calls.classes(t)
        call = t["p"] >= knob
        t["error"] = np.select([call & (t["truth"] == 0), ~call & (t["truth"] == 1), call], ["FP", "FN", "TP"], "TN")
        cols = [c for c in SPECTRUM if c in t.columns]
        print(f"\n# {rt} (knob {knob}; {len(t)} rows)\n", flush=True)
        for n in (5, 20):
            by_class(t, cols, n)
            aucs(t, cols, n)
            print()
        if "third_position_share" in t:
            insilico(t)
        print("\n## Per-taxon divergence rescaled by the gene factors (fragment tables)\n")
        rescaled(opts.records, rt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
