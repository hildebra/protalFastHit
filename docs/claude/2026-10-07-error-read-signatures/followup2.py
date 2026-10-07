#!/usr/bin/env python3
"""Second follow-up: whether the per-gene divergence's spread (cv) tells FP from FN taxa beyond their mean divergence;
the source genomes of the FN taxa's reads; the FP fragments from beyond the genus in long reads (does the rest of the
read disagree?); and the model's existing features on the same FP and FN taxa against TP and TN rows.

    python3 followup2.py --records ~/v15/err_analysis --build local/v15 --error-records ~/v15/v15_error_records
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

import fragments
import signatures

EXISTING = ["identity", "top_identity", "excess_scaled_median", "excess_high_share", "lu_gene_rate3", "lu_per_kb",
            "failed_candidate_rate", "congener_fit_share", "mean_mapq", "low_mapq_share", "fragments",
            "cluster_ani_radius", "cluster_min_ani", "cluster_mean_ani", "cluster_genomes_log10", "rep_completeness",
            "rep_contamination", "relative_spill", "genus_share", "sample_log_fragments"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--error-records", required=True, help="the unpacked v15_error_records (its samples)")
    p.add_argument("--types", default="pe,se,pb,ont")
    return p.parse_args(argv)


def per_taxon(f):
    c = f[f["counted"] & f["role"].isin(["FP", "FN"]) & f["gene"].notna()].copy()
    c["div"] = 1 - c["identity"]
    rows = []
    for (sample, taxon, role), g in c.groupby(["sample", "taxon", "role"]):
        if len(g) < 8 or g["gene"].nunique() < 4:
            continue
        gm = g.groupby("gene")["div"].mean()
        rows.append({"sample": sample, "taxon": taxon, "role": role, "div": g["div"].mean(),
                     "gene_sd": gm.std(), "cv": gm.std() / max(1e-9, gm.mean()),
                     "zero_genes": (gm == 0).mean(), "genes": len(gm), "fragments": len(g)})
    return pd.DataFrame(rows)


def spread(d):
    d = d.copy()
    d["truth"] = (d["role"] == "FN").astype(int)
    d["bin"] = pd.qcut(d["div"], 4, labels=False, duplicates="drop")
    within = [signatures.auc(g["truth"], g["cv"]) for _, g in d.groupby("bin")]
    zero = [signatures.auc(g["truth"], g["zero_genes"]) for _, g in d.groupby("bin")]
    print(f"{len(d)} taxa with >= 8 counted fragments on >= 4 genes; AUC (FN high) within quartiles of mean divergence: "
          f"cv {np.round(within, 3).tolist()}, share of genes without a difference {np.round(zero, 3).tolist()}")
    y, groups = d["truth"].to_numpy(), d["sample"].to_numpy()
    k = min(5, len(set(groups)))
    for cols in (["div"], ["div", "fragments"], ["div", "fragments", "cv"], ["div", "fragments", "cv", "zero_genes"]):
        X = np.column_stack([np.log(d[c].clip(lower=1e-4)) if c in ("div", "fragments") else d[c] for c in cols])
        s = np.zeros(len(d))
        for tr, te in GroupKFold(n_splits=k).split(X, y, groups):
            s[te] = LogisticRegression(max_iter=1000).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        print(f"  logistic, sample-grouped CV, {'+'.join(cols)}: AUC {roc_auc_score(y, s):.3f}")


def sources(f):
    for group in ("own (FN)", "wrong (FP)"):
        g = f[f["counted"] & (f["group"] == group)]
        per = g.groupby(["sample", "taxon"])["source_genome"].agg(lambda s: s.value_counts().index[0])
        print(f"  {group}: fragments by source genome {g['source_genome'].value_counts(normalize=True).round(3).to_dict()}; "
              f"taxa by their main one {per.value_counts(normalize=True).round(3).to_dict()}")


def beyond_genus(f):
    w = f[f["counted"] & (f["group"] == "wrong (FP)")]
    own = f[f["counted"] & f["group"].isin(["own (FN)", "own (TP)"])]
    cols = ["share_on_taxon", "other_taxa_records", "supplementary", "zr", "identity"]
    rows = [dict(group="FP, same genus", n=int((w["relation"] == "genus").sum()),
                 **w.loc[w["relation"] == "genus", cols].mean().round(3).to_dict()),
            dict(group="FP, beyond the genus", n=int(w["relation"].isin(["family", "order", "other"]).sum()),
                 **w.loc[w["relation"].isin(["family", "order", "other"]), cols].mean().round(3).to_dict()),
            dict(group="own reads", n=len(own), **own[cols].mean().round(3).to_dict())]
    signatures.show(pd.DataFrame(rows), index=False)
    far = w[w["relation"].isin(["family", "order", "other"])]
    if len(far):
        per = far.groupby(["sample", "taxon"]).size()
        print(f"  beyond-genus FP fragments sit on {len(per)} taxa, {per.median():.0f} per taxon (median); "
              f"the read's best share on the FP taxon < 0.5: {(far['share_on_taxon'] < 0.5).mean():.3f}")


def existing(rt, build, samples, taxa):
    t = fragments.model_tables(build, rt, lambda c: c in (["meta_sample", "taxon", "truth"] + EXISTING))
    t = t[t["meta_sample"].isin(samples)].copy()
    t["taxon"] = t["taxon"].astype(str)
    err = dict(zip(zip(taxa["sample"], taxa["taxid"].astype(str)), taxa["error"]))
    t["class"] = [err.get((s, x), "TP" if tr == 1 else "TN") for s, x, tr in zip(t["meta_sample"], t["taxon"],
                                                                                   t["truth"])]
    print(f"  {len(t)} rows of {t['meta_sample'].nunique()} samples: {t['class'].value_counts().to_dict()}")
    rows = []
    for col in EXISTING:
        if col not in t:
            continue
        x = t[col]
        r = {"feature": col}
        for cls in ("TP", "FP", "FN", "TN"):
            r[cls] = x[t["class"] == cls].median()
        a, b = t["class"].isin(["FP", "FN"]), t["class"].isin(["TP", "FP"])
        r["AUC FN vs FP"] = signatures.auc((t.loc[a, "class"] == "FN").astype(int), x[a])
        r["AUC TP vs FP"] = signatures.auc((t.loc[b, "class"] == "TP").astype(int), x[b])
        rows.append(r)
    out = pd.DataFrame(rows)
    out["strength"] = (out["AUC FN vs FP"] - 0.5).abs()
    signatures.show(out.sort_values("strength", ascending=False).drop(columns="strength").round(4), index=False)


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    for rt in opts.types.split(","):
        f = pd.read_pickle(os.path.join(opts.records, f"{rt}.fragments.pkl.gz"))
        f["group"] = signatures.group_of(f)
        taxa = fragments.load(opts.records, rt, "taxa")
        samples = {os.path.basename(os.path.dirname(os.path.dirname(p))) + ":" + os.path.basename(p)[:-len(".taxa.tsv")]
                   for p in glob.glob(os.path.join(opts.error_records, rt, "*", "*", "*.taxa.tsv"))}
        print(f"\n# {rt}\n\n## The spread of divergence over genes, beyond its mean\n", flush=True)
        spread(per_taxon(f))
        print("\n## Source genomes\n")
        sources(f)
        print("\n## FP fragments from beyond the genus\n")
        beyond_genus(f)
        print("\n## The model's existing features: FP and FN against TP and TN (medians; AUC: FN or TP high)\n")
        existing(rt, opts.build, samples, taxa)
    return 0


if __name__ == "__main__":
    sys.exit(main())
