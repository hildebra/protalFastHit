#!/usr/bin/env python3
"""Signatures of the reads behind false positives and false negatives, from fragments.py's tables.

Fragments counted for a taxon (best record at MAPQ >= 4) fall into:
  wrong (FP)     counted for an FP taxon: another species' read, by construction
  own (FN)       counted for an FN taxon, from that species: the true reads of a missed taxon
  own (TP)       counted for another taxon of the sample, from that species: a present, called taxon's true read
                 (only those reads that also touched an error taxon: a record on it, seeds on it, or its source)
  wrong (other)  counted for another taxon, from another species: a misassigned read that made no FP
Prints, per read type: where the FP fragments come from; where the FN taxa's reads went; each fragment feature's
AUC for wrong against own reads; a sample-grouped boosted model of the same; the taxa: FP against FN within bins of
counted fragments, and how the new signatures relate to the model's existing features; the genes.

    python3 signatures.py --records ~/v15/err_analysis --build local/v15 > signatures.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

import fragments

FEATURES = ["identity", "mismatches", "indels", "gap_opens", "longest_match", "clipped", "clip_bases", "at_edge",
            "rel_pos", "zu", "zt", "zr", "mapq", "mates", "mates_split", "mate_unmapped", "proper", "two_genes",
            "abs_tlen", "za_n", "za_min", "za_negative", "za_tie", "za_fn", "za_fp", "zf_n", "zf_fn", "zf_fp",
            "zf_self", "secondary", "supplementary", "other_taxa_records", "share_on_taxon", "aligned", "read_len",
            "gc", "top_trinucleotide", "homopolymer", "gene_within", "gene_between", "gene_identical",
            "gene_near_identical"]
EXISTING = ["fragments", "identity", "top_identity", "low_identity_share", "mean_mapq", "low_mapq_share", "lu_per_kb",
            "lsu_per_kb", "lu_gene_rate3", "congener_fit_share", "other_genus_fit_share", "failed_candidate_rate",
            "mate_lost_share", "excess_scaled_median", "gene_presence_ratio", "hit_gene_fraction",
            "adjacent_unlikely_share", "relative_spill", "em_own_share"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--build", required=True, help="the build's local folder (training/ and test/ tables)")
    p.add_argument("--types", default="pe,se,pb,ont")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


def auc(y, x):
    x = pd.to_numeric(pd.Series(x), errors="coerce").astype(float).to_numpy()
    y = np.asarray(y)
    ok = ~np.isnan(x)
    if ok.sum() < 10 or len(set(y[ok])) < 2:
        return np.nan
    return roc_auc_score(y[ok], x[ok])


def show(df, **kw):
    print(df.to_string(**kw), flush=True)


def group_of(f):
    own = f["relation"] == "own"
    return np.select([f["role"] == "FP", (f["role"] == "FN") & own, (f["role"] == "other") & own,
                      (f["role"] == "other") & ~own & (f["relation"] != "unknown")],
                     ["wrong (FP)", "own (FN)", "own (TP)", "wrong (other)"], default="")


def fp_sources(f, taxa):
    fn = set(zip(taxa.loc[taxa["error"] == "FN", "sample"], taxa.loc[taxa["error"] == "FN", "taxid"].astype(str)))
    w = f[f["group"] == "wrong (FP)"].copy()
    w["source_is"] = np.select(
        [w["source"] == "host", w["source_held_out"], [(s, t) in fn for s, t in zip(w["sample"], w["source_taxid"])],
         w["source_taxid"] != ""],
        ["host", "held out (novel)", "FN of the sample", "in the database (TP or unseen)"], default="unknown")
    print("Counted FP fragments by the source's relation to the FP taxon (share of fragments):")
    show(pd.crosstab(w["scenario"], w["relation"], normalize="index", margins=True).round(3))
    print("\n... by what the source species is:")
    show(pd.crosstab(w["scenario"], w["source_is"], normalize="index", margins=True).round(3))
    print("\n... by the source genome:")
    show(pd.crosstab(w["scenario"], w["source_genome"], normalize="index", margins=True).round(3))
    per = w.groupby(["sample", "taxon"]).agg(fragments=("qname", "size"), sources=("source", "nunique"),
                                             main=("relation", lambda s: s.value_counts().index[0]),
                                             novel=("source_held_out", "mean"))
    print(f"\nFP taxa with counted fragments here: {len(per)}; fragments per FP taxon: "
          f"{per['fragments'].describe(percentiles=[.5, .9]).round(1).to_dict()}")
    print(f"one source species: {(per['sources'] == 1).mean():.3f}; all fragments from held-out species: "
          f"{(per['novel'] == 1).mean():.3f}; main relation: {per['main'].value_counts(normalize=True).round(3).to_dict()}")


def fn_fate(f, taxa):
    t = taxa[taxa["error"] == "FN"].copy()
    for c in ["own_fragments", "own_best_on_taxon", "own_best_on_taxon_mapq4", "own_best_elsewhere", "own_unaligned",
              "seeded_not_aligned", "fragments_on_taxon", "best_on_taxon_mapq4"]:
        t[c] = pd.to_numeric(t[c], errors="coerce").fillna(0)
    t["counted"] = t["best_on_taxon_mapq4"]
    t["bin"] = pd.cut(t["counted"], [-1, 0, 1, 2, 5, 1e9], labels=["0", "1", "2", "3-5", ">5"])
    print("FN taxa by the fragments counted for them (best on the taxon, MAPQ >= 4):")
    show(t.groupby("scenario")["bin"].value_counts(normalize=True).unstack().round(3))
    own = t["own_fragments"].sum()
    if own:
        print(f"\nThe FN species' own fragments with a record: {int(own)}; best on itself, MAPQ >= 4: "
              f"{t['own_best_on_taxon_mapq4'].sum() / own:.3f}; best on itself below MAPQ 4: "
              f"{(t['own_best_on_taxon'] - t['own_best_on_taxon_mapq4']).sum() / own:.3f}; best elsewhere: "
              f"{t['own_best_elsewhere'].sum() / own:.3f}; no alignment: {t['own_unaligned'].sum() / own:.3f}")
    fp = set(zip(taxa.loc[taxa["error"] == "FP", "sample"], taxa.loc[taxa["error"] == "FP", "taxid"].astype(str)))
    fnset = set(zip(t["sample"], t["taxid"].astype(str)))
    s = f[[(a, b) in fnset for a, b in zip(f["sample"], f["source_taxid"])]]
    if len(s):
        s = s.assign(went=np.select([s["taxon"] == "", s["relation"] == "own",
                                     [(a, b) in fp for a, b in zip(s["sample"], s["taxon"])]],
                                    ["unaligned", "itself", "an FP taxon"], default="another taxon"))
        print("\nThe FN species' fragments in the SAMs (capped per reason in newer builds), where their best record is, "
              "by the relation of that taxon:")
        show(pd.crosstab(s["went"], s["relation"]).astype(int))
        swaps = s[s["went"] == "an FP taxon"].groupby(["sample", "taxon"]).size()
        print(f"FP taxa fed by an FN species of the same sample: {len(swaps)} of {len(fp)} FP taxa")
        cnt = s[(s["went"] == "itself")]
        if len(cnt):
            print(f"FN own fragments on itself: counted (MAPQ >= 4) {cnt['counted'].mean():.3f}; of those below MAPQ 4, "
                  f"ZA negative (another taxon fits better) {cnt.loc[~cnt['counted'], 'za_negative'].mean():.3f}, "
                  f"ZA tie {cnt.loc[~cnt['counted'], 'za_tie'].mean():.3f}")


def contrasts(f, seed):
    c = f[f["counted"] & (f["group"] != "")]
    print(pd.Series(c["group"]).value_counts().to_string())
    w = c["group"].str.startswith("wrong (FP)")
    rows = []
    for col in FEATURES:
        if col not in c or c[col].isna().all():
            continue
        x = c[col].astype(float)
        row = {"feature": col}
        for g in ["wrong (FP)", "wrong (other)", "own (FN)", "own (TP)"]:
            row[g] = x[c["group"] == g].mean()
        fp_or_own = c["group"].isin(["wrong (FP)", "own (FN)", "own (TP)"])
        row["AUC FP vs own"] = auc(w[fp_or_own], x[fp_or_own])
        sub = c["group"].isin(["wrong (FP)", "own (FN)"])
        row["AUC FP vs own FN"] = auc(w[sub], x[sub])
        sub = c["group"].isin(["wrong (FP)", "own (TP)"])
        row["AUC FP vs own TP"] = auc(w[sub], x[sub])
        rows.append(row)
    out = pd.DataFrame(rows)
    out["strength"] = (out["AUC FP vs own"] - 0.5).abs()
    show(out.sort_values("strength", ascending=False).drop(columns="strength").round(3), index=False)

    d = c[c["group"].isin(["wrong (FP)", "own (FN)", "own (TP)"])]
    y = (d["group"] == "wrong (FP)").astype(int).to_numpy()
    cols = [x for x in FEATURES if x in d and not d[x].isna().all()]
    groups = d["sample"].to_numpy()
    k = min(5, len(set(groups)))
    if k >= 2 and y.sum() >= 20 and (1 - y).sum() >= 20:
        s = np.zeros(len(d))
        for tr, te in GroupKFold(n_splits=k).split(d, y, groups):
            m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, random_state=seed)
            m.fit(d[cols].astype(float).iloc[tr], y[tr])
            s[te] = m.predict_proba(d[cols].astype(float).iloc[te])[:, 1]
        ident = auc(y, -d["identity"].astype(float))
        print(f"\nboosted model, wrong (FP) against own reads, sample-grouped {k}-fold: AUC {roc_auc_score(y, s):.3f} "
              f"(identity alone {ident:.3f}); {len(d)} fragments")
        for name, drop in (("without identity, mismatches, longest match", ["identity", "mismatches",
                                                                              "longest_match"]),
                           ("without the ZA/ZF tags", [x for x in cols if x.startswith(("za_", "zf_"))])):
            cc = [x for x in cols if x not in drop]
            s2 = np.zeros(len(d))
            for tr, te in GroupKFold(n_splits=k).split(d, y, groups):
                m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, random_state=seed)
                m.fit(d[cc].astype(float).iloc[tr], y[tr])
                s2[te] = m.predict_proba(d[cc].astype(float).iloc[te])[:, 1]
            print(f"  {name}: AUC {roc_auc_score(y, s2):.3f}")


def taxon_level(f, taxa, table, seed):
    c = f[f["counted"] & f["role"].isin(["FP", "FN"])].copy()
    flags = ["za_negative", "za_tie", "za_fn", "zf_fn", "zf_self", "mates_split", "mate_unmapped", "clipped",
             "at_edge", "two_genes"]
    for x in flags:
        c[x] = c[x].astype(float)
    agg = c.groupby(["sample", "taxon", "role"]).agg(
        n=("qname", "size"), identity=("identity", "mean"), zu=("zu", "mean"), zt=("zt", "mean"),
        zu0=("zu", lambda s: (s == 0).mean()), gene_identical=("gene_identical", "mean"),
        genes=("gene", "nunique"), share_on_taxon=("share_on_taxon", "mean"),
        **{x: (x, "mean") for x in flags}).reset_index()
    agg["truth"] = (agg["role"] == "FN").astype(int)
    agg["bin"] = pd.cut(agg["n"], [0, 1, 2, 5, 1e9], labels=["1", "2", "3-5", ">5"])
    print(pd.crosstab(agg["role"], agg["bin"]).to_string())
    cols = ["identity", "zu", "zt", "zu0", "gene_identical", "genes", "share_on_taxon"] + flags
    rows = []
    for x in cols:
        row = {"feature": x, "FP mean": agg.loc[agg["truth"] == 0, x].mean(),
               "FN mean": agg.loc[agg["truth"] == 1, x].mean(), "AUC (FN high)": auc(agg["truth"], agg[x])}
        for b, g in agg.groupby("bin", observed=True):
            row[f"AUC n={b}"] = auc(g["truth"], g[x])
        rows.append(row)
    out = pd.DataFrame(rows)
    show(out.round(3), index=False)
    if table is not None:
        t = table.rename(columns={"meta_sample": "sample"})
        t["taxon"] = t["taxon"].astype(str)
        j = agg.merge(t, on=["sample", "taxon"], how="left", suffixes=("", "_model"))
        print(f"\njoined to the model's table: {j['fragments'].notna().sum()} of {len(j)} taxa; Spearman correlation of "
              "the new signatures with the model's features (FP and FN taxa):")
        ex = [e if e in j else e + "_model" for e in EXISTING]
        ex = [e for e in ex if e in j]
        corr = j[cols + ex].astype(float).corr(method="spearman").loc[cols, ex]
        show(corr.round(2))


def genes(f):
    c = f[f["counted"] & f["group"].isin(["wrong (FP)", "own (FN)", "own (TP)"])]
    by = pd.crosstab(c["gene"], c["group"])
    if "wrong (FP)" not in by or by.empty:
        return
    own = by.drop(columns="wrong (FP)").sum(axis=1)
    by["FP share"] = by["wrong (FP)"] / by["wrong (FP)"].sum()
    by["own share"] = own / own.sum()
    by["enrichment"] = (by["FP share"] / by["own share"].replace(0, np.nan)).round(2)
    gi = c.groupby("gene")["gene_identical"].first()
    by["identical_share"] = gi
    by = by.sort_values("enrichment", ascending=False)
    top = by["FP share"].sort_values(ascending=False)
    print(f"genes with FP fragments: {(by['wrong (FP)'] > 0).sum()}; the 10 most FP-rich hold "
          f"{top.head(10).sum():.3f} of the FP fragments (own reads: {by.loc[top.head(10).index, 'own share'].sum():.3f})")
    print(f"Spearman(enrichment, congeners' identical share): "
          f"{by['enrichment'].corr(by['identical_share'], method='spearman'):.2f}")
    show(by.head(12).round(4))


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    for rt in opts.types.split(","):
        f = pd.read_pickle(os.path.join(opts.records, f"{rt}.fragments.pkl.gz"))
        taxa = fragments.load(opts.records, rt, "taxa")
        taxa["taxid"] = taxa["taxid"].astype(str)
        f["group"] = group_of(f)
        table = fragments.model_tables(opts.build, rt, lambda c: c in (["meta_sample", "taxon", "truth"] + EXISTING))
        print(f"\n# {rt}: {taxa['sample'].nunique()} samples, {int((taxa['error'] == 'FP').sum())} FP and "
              f"{int((taxa['error'] == 'FN').sum())} FN taxa; {len(f)} fragments\n", flush=True)
        print("## Where the FP taxa's counted fragments come from\n")
        fp_sources(f, taxa)
        print("\n## Where the FN taxa's reads went\n")
        fn_fate(f, taxa)
        print("\n## Read level: counted fragments, wrong (FP) against own reads (means; AUC > 0.5: higher in FP reads)\n")
        contrasts(f, opts.seed)
        print("\n## Taxon level: FP against FN taxa, from their counted fragments (AUC > 0.5: higher in FN)\n")
        taxon_level(f, taxa, table, opts.seed)
        print("\n## Genes: where the FP fragments lie against the own reads\n")
        genes(f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
