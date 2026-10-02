#!/usr/bin/env python3
"""The present taxa of the r226 v3 long-read training data that the model misses with species held out (p_species
below 0.5) although they have more than 100 fragments, next to the found ones with as many: their features, from the
training tables. Usage: strains_missed.py V3_DIR"""
import os
import sys

import pandas as pd

SHOWN = ["fragments", "identity", "top_identity", "hit_gene_fraction", "excess_median", "excess_high_share",
         "conserved_fast_depth_ratio", "conserved_hit_share", "congener_fit_share", "low_mapq_share"]


def main():
    v3 = sys.argv[1]
    pd.set_option("display.width", 250)
    for t in ("pb", "ont"):
        pred = pd.read_csv(os.path.join(v3, f"trained_model_{t}.predictions.tsv.gz"), sep="\t",
                           usecols=["meta_sample", "taxon", "taxon_name", "truth", "meta_rep_genome", "p_species"])
        table = pd.read_csv(os.path.join(v3, "training", f"training_data_{t}.tsv"), sep="\t",
                            usecols=["meta_sample", "taxon"] + SHOWN)
        df = pred.merge(table, on=["meta_sample", "taxon"])
        deep = df[(df["truth"] == 1) & (df["fragments"] > 100)]
        missed = deep[deep["p_species"] < 0.5]
        print(f"## {t}: {len(missed)} of {len(deep)} present taxa with more than 100 fragments missed "
              f"({int((missed['meta_rep_genome'] == 0).sum())} strains, i.e. not the representative)")
        print(missed.sort_values("fragments", ascending=False)[["meta_sample", "taxon_name", "meta_rep_genome", "p_species"]
                                                              + SHOWN].round(4).to_string(index=False))
        print("medians:")
        groups = {"missed": missed, "found strains": deep[(deep["p_species"] >= 0.5) & (deep["meta_rep_genome"] == 0)],
                  "found representatives": deep[(deep["p_species"] >= 0.5) & (deep["meta_rep_genome"] == 1)]}
        print(pd.DataFrame({name: g[SHOWN].median() for name, g in groups.items()}).T.round(4).to_string())
        print()


if __name__ == "__main__":
    main()
