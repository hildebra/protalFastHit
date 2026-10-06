#!/usr/bin/env python3
"""The highest-scored false positives of the soil samples, with the taxid that names their genes in the SAM files
(references <taxid>_<geneid>) and where the build keeps each sample (training/ or test/ of its collection folder).

    python3 top_fp_list.py --rows v13_rows.tsv.gz --out top_false_positives.tsv
"""
import argparse

import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--rows", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--top", type=int, default=30)
opts = ap.parse_args()
df = pd.read_csv(opts.rows, sep="\t", low_memory=False)
df = df[df.scenario.isin(["soil", "soil_shallow"])]
fp = df[df.outcome == "FP"].copy()
fp["collection"] = fp.set.map({"hold-in species held out": "training", "hold-out": "test"})
fp["sample_key"] = fp.collection + "/" + fp.meta_sample
fp["fp_samples"] = fp.groupby(["read_type", "scenario", "taxon"]).sample_key.transform("nunique")
fp["point"] = fp.meta_sample.str.replace(r"_s_\d+$", "", regex=True)
cols = ["read_type", "scenario", "collection", "point", "meta_sample", "taxon", "taxon_name", "p", "fragments",
        "identity", "excess_scaled_median", "low_mapq_share", "congener_fit_share", "fp_samples"]
top = (fp.sort_values("p", ascending=False).groupby(["read_type", "scenario"], group_keys=False)
       .head(opts.top)[cols].sort_values(["read_type", "scenario", "p"], ascending=[True, True, False]))
top.to_csv(opts.out, sep="\t", index=False, float_format="%.4g")
print(top[top.read_type == "pe"].head(12).to_string(index=False))
