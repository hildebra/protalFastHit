"""A taxon beside a held-out congener (meta_novel_distance <= 0.15) with 3 or more fragments: present (its own reads
plus the congener's) or absent (only the congener's). How well each feature, and the model, tells the two apart
(AUC), and the features' medians in four groups: present and called (TP), present and missed with most of its reads
the congener's (swamped FN, error_reads), present and missed otherwise, absent and called (FP). Also the missed
taxa's own reads that the error-reads table counts (own_best_on_taxon): what an oracle count of the taxon's own
reads would have to work with.

usage: twin_separability.py <build dir> <read type>
"""
import sys
import numpy as np
import pandas as pd

D, T = sys.argv[1], sys.argv[2]
SUF = "" if T == "pe" else "_" + T
pd.set_option("display.width", 250)
FEAT = ["fragments", "identity", "top_identity", "low_identity_share", "lu_per_kb", "excess_median",
        "excess_scaled_median", "excess_high_share", "congener_fit_share", "unexpected_congener_fit_share",
        "conserved_fast_kept_ratio", "gene_presence_ratio", "hit_gene_fraction", "allele_explained_share",
        "low_mapq_share", "relative_close_share", "read_consensus_share"]


def auc(pos, neg):
    """P(a present taxon's value > an absent one's), ties half."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    pos, neg = pos[~np.isnan(pos)], neg[~np.isnan(neg)]
    if not len(pos) or not len(neg):
        return np.nan
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank().to_numpy()
    return (ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


calls = pd.read_csv(f"{D}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t",
                    usecols=["meta_sample", "meta_novel_distance", "taxon_name", "truth", "set", "p", "call"])
calls = calls[calls["set"] == "training"]
hdr = pd.read_csv(f"{D}/work/training/training_data{SUF}.tsv", sep="\t", nrows=0).columns
t = pd.read_csv(f"{D}/work/training/training_data{SUF}.tsv", sep="\t",
                usecols=["meta_sample", "taxon_name"] + [f for f in FEAT if f in hdr]).drop_duplicates(
    ["meta_sample", "taxon_name"])
c = calls.merge(t, on=["meta_sample", "taxon_name"], how="left")
dist = pd.to_numeric(c["meta_novel_distance"], errors="coerce")
c = c[(dist <= 0.15) & (c["fragments"] >= 3)].copy()
err = pd.read_csv(f"{D}/model_logs/error_reads/{T}/taxa.tsv.gz", sep="\t")
err = err[(err["set"] == "training") & (err["error"] == "FN")].copy()


def own_share(r):
    tot = own = 0
    for p in str(r["fragments_on_taxon_from"]).split(","):
        if ":" in p:
            name, n = p.rsplit(":", 1)
            tot += int(n)
            own += int(n) if name == r["taxon_name"] else 0
    return own / tot if tot else np.nan


err["own_share"] = err.apply(own_share, axis=1)
c = c.merge(err[["sample", "taxon_name", "own_share", "own_best_on_taxon"]], left_on=["meta_sample", "taxon_name"],
            right_on=["sample", "taxon_name"], how="left")
present, called = c["truth"] == 1, c["call"] == 1
c["group"] = np.select([present & called, present & ~called & (c["own_share"] < 0.5), present & ~called,
                        ~present & called, ~present & ~called],
                       ["TP", "FN swamped", "FN other", "FP", "TN"], default="")
print(f"######## {T}: taxa beside a held-out congener within 0.15 with 3 or more fragments (training, species held out)")
print(c["group"].value_counts().to_string())
print("\n== AUC present vs absent (all of them), and present-but-missed (FN) vs absent-and-called (FP)")
rows = []
fnm, fpm = c["group"].str.startswith("FN"), c["group"] == "FP"
for f in ["p"] + [f for f in FEAT if f in c.columns]:
    rows.append({"feature": f, "AUC present/absent": auc(c.loc[present, f], c.loc[~present, f]),
                 "AUC FN/FP": auc(c.loc[fnm, f], c.loc[fpm, f])})
print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print("\n== medians by group")
print(c.groupby("group")[[f for f in FEAT if f in c.columns] + ["p"]].median().T.to_string(
    float_format=lambda v: f"{v:.4f}"))
sw = c[c["group"] == "FN swamped"]
print(f"\n== swamped FN: own reads on the taxon (own_best_on_taxon) median {sw['own_best_on_taxon'].median():.0f}, "
      f"10% {sw['own_best_on_taxon'].quantile(0.1):.0f}, 90% {sw['own_best_on_taxon'].quantile(0.9):.0f}; own share median "
      f"{sw['own_share'].median():.3f}; with 10 or more own reads: {int((sw['own_best_on_taxon'] >= 10).sum())} of {len(sw)}")
