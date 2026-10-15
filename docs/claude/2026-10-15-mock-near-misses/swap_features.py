#!/usr/bin/env python3
"""What a taxon's probability owes to each feature: its row of one protal training dump (<profile>.truth_annotated)
scored by the PMML model (scripts/model_pmml.py, as protal scores it), then again with one feature at a time, and with
named groups of features, set to a reference: the same taxon's row in another dump (--other DUMP), or the median of
the training table's present rows with fragments in [--min-fragments, --max-fragments] (--training TABLE). Prints the
features sorted by the change in probability. Needs numpy and pandas.

usage: swap_features.py MODEL.xml DUMP TAXID (--other DUMP | --training TABLE) [--top 20]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import model_pmml  # noqa: E402

GROUPS = {
    "sample": ["sample_log_fragments", "sample_log_taxa", "sample_low_identity", "sample_identity",
               "sample_read_length", "sample_read_length_sd"],
    "read_length": ["top_identity_bound", "allele_top_identity_bound", "aligned_length_ratio", "sample_read_length",
                    "sample_read_length_sd"],
    "coverage": ["depth_cv", "breadth_ratio", "gene_presence_ratio", "hit_gene_fraction", "gene_dispersion",
                 "gene_divergence_dispersion"],
    "identity": ["identity", "top_identity", "allele_identity", "allele_top_identity", "low_identity_share",
                 "excess_median", "excess_scaled_median", "excess_high_share", "conserved_mismatch_ratio",
                 "conserved_mismatch_rate", "nonsynonymous_conserved_rate", "third_position_share",
                 "nonsynonymous_share", "fixed_difference_rate"],
    "alleles": ["allele_explained_share", "allele_identity_gain", "polymorphic_known_share", "polymorphic_site_rate",
                "ancestry_fixed_gain", "ancestry_fixed_agreement"],
    "congeners": ["congener_fit_share", "unexpected_congener_fit_share", "low_mapq_share", "mean_mapq", "uniqueness",
                  "seed_crowding", "other_genus_fit_share"],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model")
    ap.add_argument("dump")
    ap.add_argument("taxid")
    ap.add_argument("--other")
    ap.add_argument("--training")
    ap.add_argument("--min-fragments", type=float, default=200)
    ap.add_argument("--max-fragments", type=float, default=20000)
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()
    if bool(args.other) == bool(args.training):
        sys.exit("give one of --other and --training")

    model = model_pmml.load_model(args.model)
    features = model.features
    table = pd.read_csv(args.dump, sep="\t", low_memory=False)
    row = table[table["taxon"].astype(str) == args.taxid]
    if row.empty:
        sys.exit(f"taxon {args.taxid} not in {args.dump}")
    x = row.iloc[[0]][features].astype(float).reset_index(drop=True)
    if args.other:
        other = pd.read_csv(args.other, sep="\t", low_memory=False)
        orow = other[other["taxon"].astype(str) == args.taxid]
        if orow.empty:
            sys.exit(f"taxon {args.taxid} not in {args.other}")
        reference = orow.iloc[0][features].astype(float)
        what = f"its row in {os.path.basename(args.other)} (P {float(orow.iloc[0]['probability']):.3f})"
    else:
        train = pd.read_csv(args.training, sep="\t", low_memory=False)
        present = train[(train["truth"] == 1) & train["fragments"].between(args.min_fragments, args.max_fragments)]
        reference = present[features].astype(float).median()
        what = f"the median of {len(present)} present training rows"
    p0 = float(model.predict(x)[0])
    print(f"# {os.path.basename(args.dump)} taxon {args.taxid} {row.iloc[0]['taxon_name']}: protal's P "
          f"{float(row.iloc[0]['probability']):.4f}, the model's here {p0:.4f}; reference: {what}")

    def score(columns):
        trial = x.copy()
        for c in columns:
            trial[c] = reference[c]
        return float(model.predict(trial)[0])

    print("group\tP\tdP")
    for name, columns in GROUPS.items():
        columns = [c for c in columns if c in features]
        if columns:
            p = score(columns)
            print(f"{name} ({len(columns)})\t{p:.4f}\t{p - p0:+.4f}")
    p_all = score(features)
    print(f"all\t{p_all:.4f}\t{p_all - p0:+.4f}")

    rows = []
    for f in features:
        if np.isclose(float(x.at[0, f]), float(reference[f]), equal_nan=True):
            continue
        p = score([f])
        rows.append((f, float(x.at[0, f]), float(reference[f]), p - p0))
    rows.sort(key=lambda r: -abs(r[3]))
    print("feature\tvalue\treference\tdP")
    for f, own, ref, dp in rows[:args.top]:
        print(f"{f}\t{own:.6g}\t{ref:.6g}\t{dp:+.4f}")


if __name__ == "__main__":
    sys.exit(main())
