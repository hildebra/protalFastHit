#!/usr/bin/env python3
"""Checks on the FN taxa behind the signatures: their counted and failed reads per taxon, whose reads seed on them and
fail, the FN taxa with many counted fragments, and the most recurrent of those in the model's table (pe).

    python3 checks.py --records ~/v15/err_analysis --build local/v15 > checks.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

import fragments

TABLE_COLUMNS = ["meta_sample", "taxon_name", "truth", "fragments", "identity", "top_identity", "excess_scaled_median",
                 "lu_per_kb", "lu_gene_rate3", "hit_gene_fraction", "rep_completeness", "cluster_genomes_log10",
                 "cluster_min_ani", "failed_candidate_rate", "low_mapq_share"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--types", default="pe,se,pb,ont")
    return p.parse_args(argv)


def numeric(t, cols):
    for c in cols:
        t[c] = pd.to_numeric(t[c], errors="coerce").fillna(0)
    return t


def fn_counts(rt, t):
    fn = numeric(t[t["error"] == "FN"].copy(), ["own_best_on_taxon", "own_best_on_taxon_mapq4", "own_best_elsewhere",
                                                "seeded_not_aligned", "best_on_taxon_mapq4"])
    print(f"{rt}: {len(fn)} FN taxa; counted median {fn['best_on_taxon_mapq4'].median():.0f}, seeded on it and not "
          f"aligned median {fn['seeded_not_aligned'].median():.0f} (more than counted: "
          f"{(fn['seeded_not_aligned'] > fn['best_on_taxon_mapq4']).mean():.3f}); own reads best elsewhere: "
          f"{(fn['own_best_elsewhere'] > 0).mean():.3f}; own reads on it below MAPQ 4: "
          f"{((fn['own_best_on_taxon'] - fn['own_best_on_taxon_mapq4']) > 0).mean():.3f}")
    print("  medians by scenario:", fn.groupby("scenario")[["best_on_taxon_mapq4", "seeded_not_aligned",
                                                            "own_best_elsewhere"]].median().to_dict("index"))


def seeded(rt, f):
    s = f[f["xe"].fillna("").str.contains("seeded:")].copy()
    s["seeded_on"] = s["xe"].str.findall(r"seeded:(\d+)")
    s = s.explode("seeded_on").reset_index(drop=True)
    s["whose"] = np.where(s["seeded_on"] == s["source_taxid"], "the FN species' own read", "another species' read")
    s["fate"] = np.where(s["taxon"] == "", "aligned nowhere", "aligned elsewhere")
    print(f"{rt}: fragments that seeded on an FN taxon and did not align to it: {len(s)}")
    print(pd.crosstab(s["whose"], s["fate"]).to_string())
    own = s[s["whose"] == "the FN species' own read"]
    if len(own):
        per = own.groupby(["sample", "seeded_on"]).size()
        print(f"  own failed reads per FN taxon with any: median {per.median():.0f}, 90% {per.quantile(.9):.0f}; "
              f"{len(per)} FN taxa")
    others = s.loc[s["whose"] == "another species' read", "source_held_out"].mean()
    print(f"  the others' failed reads from held-out species: {others:.3f}")


def big_fn(rt, t):
    t = numeric(t.copy(), ["own_best_elsewhere", "best_on_taxon_mapq4"])
    fn = t[t["error"] == "FN"]
    big = fn[fn["best_on_taxon_mapq4"] >= 20]
    fp = set(zip(t.loc[t["error"] == "FP", "sample"], t.loc[t["error"] == "FP", "taxon_name"]))
    went = [("FP" if (s, x.split(",")[0].rsplit(":", 1)[0]) in fp else "not FP")
            for s, x in zip(big["sample"], big["own_best_elsewhere_on"].fillna("")) if x]
    print(f"{rt}: FN taxa with >= 20 counted fragments: {len(big)} of {len(fn)}; {big['scenario'].value_counts().to_dict()}; "
          f"p median {pd.to_numeric(big['p']).median():.2f}; own reads mostly elsewhere: "
          f"{(big['own_best_elsewhere'] > big['best_on_taxon_mapq4']).mean():.2f}; the taxon most of those went to: "
          f"{pd.Series(went).value_counts().to_dict()}")


def recurrent(build, records):
    t = fragments.model_tables(build, "pe", lambda c: c in TABLE_COLUMNS)
    er = numeric(fragments.load(records, "pe", "taxa"), ["best_on_taxon_mapq4"])
    big = er[(er["error"] == "FN") & (er["best_on_taxon_mapq4"] >= 20)]
    rec = big["taxon_name"].value_counts()
    print("pe FN species with >= 20 counted fragments, most recurrent:", rec.head(8).to_dict())
    for name in rec.head(5).index:
        rows = t[t["taxon_name"] == name]
        print(f"\n{name}: {len(rows)} rows, {int(rows['truth'].sum())} present; medians by truth")
        print(rows.groupby("truth")[TABLE_COLUMNS[3:]].median().round(3).T.to_string())
    cls = dict(zip(zip(er["sample"], er["taxon_name"]), er["error"]))
    t["class"] = [cls.get((s, n), "TP" if tr == 1 else "TN") for s, n, tr in zip(t["meta_sample"], t["taxon_name"],
                                                                                  t["truth"])]
    sub = t[(t["fragments"] >= 20) & t["class"].isin(["TP", "FN", "FP"])]
    print(f"\nrows with >= 20 fragments {sub['class'].value_counts().to_dict()}, medians by class:")
    print(sub.groupby("class")[TABLE_COLUMNS[3:]].median().round(3).T.to_string())


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    for rt in opts.types.split(","):
        t = fragments.load(opts.records, rt, "taxa")
        f = pd.read_pickle(os.path.join(opts.records, f"{rt}.fragments.pkl.gz"))
        print(f"\n# {rt}\n")
        fn_counts(rt, t)
        seeded(rt, f)
        big_fn(rt, t)
    print("\n# pe: the recurrent FN species with many fragments\n")
    recurrent(opts.build, opts.records)
    return 0


if __name__ == "__main__":
    sys.exit(main())
