#!/usr/bin/env python3
"""Soil rows by divergence from the reference (soil_errors.py --save rows): present taxa's miss rate and absent
taxa's false-positive rate by excess_scaled_median (the reads' divergence beyond their base qualities, scaled by gene
conservation: a genome-wide divergence estimate) and by the GTDB cluster's minimum intra-species ANI, per strain kind.
Shows where real strains and the congeners of novel species overlap.

    python3 soil_divergence.py --rows v13_rows.tsv.gz
"""
import argparse

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", required=True)
    ap.add_argument("--read-types", default="pe,pb")
    opts = ap.parse_args()
    pd.set_option("display.width", 250)
    df = pd.read_csv(opts.rows, sep="\t", low_memory=False)
    df["kind"] = np.where(df.truth == 1, df.strain, np.where(df.meta_novel_level.astype(str) == "species",
                                                             "absent: congener of a novel species", "absent: other"))
    ex = pd.cut(df.excess_scaled_median, [-1, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 1])
    ani = pd.cut(df.cluster_min_ani, [-2, 0, 96, 97, 98, 99, 101], labels=["unknown/one genome", "<96", "96-97", "97-98", "98-99", ">=99"])
    fr = pd.cut(df.fragments, [0, 3, 10, 1e9], right=False, labels=["1-2", "3-9", ">=10"])
    for rt in opts.read_types.split(","):
        for sc in ["soil", "soil_shallow"]:
            g = df[(df.read_type == rt) & (df.scenario == sc)]
            if not len(g):
                continue
            err = np.where(g.truth == 1, g.outcome == "FN", g.outcome == "FP")
            g = g.assign(err=err)
            print(f"\n=== {rt} {sc} (hold-in species held out and hold-out together), >= 3 fragments")
            h = g[g.fragments >= 3]
            t = h.pivot_table(index=ex[h.index], columns="kind", values="err", observed=True,
                              aggfunc=lambda s: f"{100 * s.mean():5.1f}% ({len(s)})")
            print("error rate (rows) by excess_scaled_median:")
            print(t.to_string())
            pres = g[g.truth == 1]
            print("present: FN rate (rows) by the cluster's minimum ANI and fragments:")
            print(pres.pivot_table(index=ani[pres.index], columns=[pres.strain, fr[pres.index]], values="err",
                                   observed=True, aggfunc=lambda s: f"{100 * s.mean():.1f} ({len(s)})").to_string())


if __name__ == "__main__":
    main()
