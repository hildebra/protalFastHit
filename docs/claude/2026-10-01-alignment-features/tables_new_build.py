#!/usr/bin/env python3
"""Training tables (collect_training_data.py's format) from the dumps of profile_new_build.sh: each dump's rows with
the meta columns of the V2 tables (per sample, and per taxon where V2's run saw the taxon too; otherwise the domain
from the training database's lineages and the per-taxon columns empty).

    python3 tables_new_build.py NEW_OUT
"""
import os
import sys

import pandas as pd

V2 = os.path.expanduser(os.environ.get("V2", "~/tune/V2"))
OUT = os.path.expanduser(sys.argv[1])
TABLE = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv", "ont": "training_data_ont.tsv"}
SAMPLE_META = ["meta_design", "meta_sample", "meta_read_length", "meta_read_pairs", "meta_novel_species", "meta_novel_levels",
               "meta_read_type"]
TAXON_META = ["meta_domain", "meta_novel_congener", "meta_rep_genome", "meta_novel_level", "meta_relative_rank",
              "meta_neighbour_rank"]
META = ["meta_design", "meta_sample", "meta_read_length", "meta_read_pairs", "meta_domain", "meta_novel_species",
        "meta_novel_congener", "meta_rep_genome", "meta_novel_levels", "meta_novel_level", "meta_relative_rank",
        "meta_neighbour_rank", "meta_read_type"]

domain_of = {}
with open(f"{V2}/training_db/genome2tiid.tsv") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4:
            lineage = f[3].split(";")
            domain_of["s__" + lineage[-1][3:]] = lineage[0][3:]

for split in ("training", "test"):
    rows = []
    with open(f"{OUT}/{split}/samples.map") as fh:
        header = None
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "#SAMPLEID":
                header = [x.lstrip("#") for x in f]
            elif not f[0].startswith("#") and header:
                rows.append(dict(zip(header, f)))
    by_type = {}
    for r in rows:
        by_type.setdefault(r["READ_TYPE"], []).append(r)
    for rt, samples in by_type.items():
        old = pd.read_csv(f"{V2}/{split}/{TABLE[rt]}", sep="\t", low_memory=False)
        per_sample = old.drop_duplicates("meta_sample").set_index("meta_sample")[[c for c in SAMPLE_META if c != "meta_sample"]]
        per_taxon = old.drop_duplicates(["meta_sample", "taxon_name"]).set_index(["meta_sample", "taxon_name"])[TAXON_META]
        frames = []
        for r in samples:
            dump = r["PROFILE"] + ".truth_annotated"
            d = pd.read_csv(dump, sep="\t", float_precision="round_trip", low_memory=False)
            sample = r["SAMPLEID"]
            meta = pd.DataFrame(index=d.index)
            meta["meta_sample"] = sample
            for c in per_sample.columns:
                meta[c] = per_sample.loc[sample, c]
            keyed = per_taxon.reindex(pd.MultiIndex.from_arrays([[sample] * len(d), d["taxon_name"]]))
            for c in TAXON_META:
                meta[c] = keyed[c].to_numpy()
            missing = meta["meta_domain"].isna()
            meta.loc[missing, "meta_domain"] = d.loc[missing, "taxon_name"].map(domain_of).fillna("unknown")
            frames.append(pd.concat([meta[META], d], axis=1))
        table = pd.concat(frames, ignore_index=True)
        os.makedirs(f"{OUT}/{split}", exist_ok=True)
        table.to_csv(f"{OUT}/{split}/{TABLE[rt]}", sep="\t", index=False)
        matched = table["meta_relative_rank"].notna().mean()
        print(f"{split} {rt}: {len(table)} taxa in {len(samples)} samples ({int(table['truth'].astype(str).str.lower().isin(['1', 'true']).sum())} "
              f"present); V2 had {len(old)}; {100 * matched:.1f}% of the rows matched to V2's", flush=True)
