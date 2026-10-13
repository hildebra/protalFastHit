"""v19 against v21: FN and FP by genome class, split by whether the species has strain alleles in v21's training
database; and how well v21's allele/polymorphic/ancestry features separate missed-strain-like rows from FP-like ones.

usage: alleles_classes.py <local dir> <threshold>
"""
import sys
import numpy as np
import pandas as pd

L = sys.argv[1]
TH = float(sys.argv[2])
V19 = f"{L}/v19/protal0.7.9_r226_v19"
V21 = f"{L}/v21/protal0.7.9_r226_v21"
SUF = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
FEATS = ["allele_explained_share", "allele_identity_gain", "polymorphic_known_share", "polymorphic_novel_share",
         "ancestry_fixed_gain", "ancestry_agreement", "ancestry_congener_share", "identity"]


def auc(pos, neg):
    """P(pos > neg) + 0.5 P(tie)."""
    pos = np.asarray(pos, float)
    neg = np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank().to_numpy()
    rp = ranks[:len(pos)].sum()
    return (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def gclass(df):
    rep = pd.to_numeric(df["meta_rep_genome"], errors="coerce") == 1
    ins = pd.to_numeric(df["meta_insilico_strain"], errors="coerce") == 1
    return np.where(rep, "rep", np.where(ins, "insilico", "real"))


for t, suf in SUF.items():
    cols = ["taxon_name", "fragments", "allele_copy_share"]
    parts = [pd.read_csv(f"{V21}/work/{d}/training_data{suf}.tsv", sep="\t", usecols=cols) for d in ("training", "test")]
    tab = pd.concat(parts)
    tab = tab[tab["fragments"] >= 3]
    has = tab.groupby("taxon_name")["allele_copy_share"].max()
    print(f"== {t}: species seen with >= 3 fragments {len(has)}, with allele copies (max share > 0.2) "
          f"{(has > 0.2).mean():.3f}")
    rows = []
    for v, base in (("v19", V19), ("v21", V21)):
        c = pd.read_csv(f"{base}/model_logs/trained_model{suf}.calls.tsv.gz", sep="\t",
                        usecols=["meta_scenario", "meta_rep_genome", "meta_insilico_strain", "meta_novel_congener",
                                 "taxon_name", "truth", "set", "p"], dtype={"meta_scenario": str})
        c["call"] = c["p"] >= TH
        c["alleles"] = c["taxon_name"].map(has).fillna(-1)
        c["agrp"] = np.where(c["alleles"] < 0, "unseen", np.where(c["alleles"] > 0.2, "alleles", "none"))
        pres = c[c["truth"] == 1].copy()
        pres["cls"] = gclass(pres)
        for (s, cl, g), d in pres.groupby(["set", "cls", "agrp"]):
            rows.append((v, s, cl, g, len(d), int((~d["call"]).sum())))
        ab = c[c["truth"] == 0].copy()
        ab["cls"] = np.where(pd.to_numeric(ab["meta_novel_congener"], errors="coerce") == 1, "abs_beside_heldout", "abs_other")
        for (s, cl, g), d in ab.groupby(["set", "cls", "agrp"]):
            rows.append((v, s, cl, g, len(d), int(d["call"].sum())))
    r = pd.DataFrame(rows, columns=["build", "set", "class", "alleles", "n", "err"])
    r["rate"] = r["err"] / r["n"]
    w = r.pivot_table(index=["set", "class", "alleles"], columns="build", values=["n", "err", "rate"], dropna=False)
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w[["n_v19", "err_v19", "rate_v19", "n_v21", "err_v21", "rate_v21"]].reset_index()
    with pd.option_context("display.width", 200, "display.max_rows", 200, "display.float_format", "{:.4f}".format):
        print(w.to_string(index=False))

    # separation in v21 (training table, species held out p from the calls): present real strains with alleles
    # against absent rows beside a held-out congener, among the rows the model is unsure of or wrong about
    tr = pd.read_csv(f"{V21}/work/training/training_data{suf}.tsv", sep="\t",
                     usecols=["meta_sample", "taxon", "truth", "meta_rep_genome", "meta_insilico_strain",
                              "meta_novel_congener", "fragments", "allele_copy_share"] + FEATS)
    pr = pd.read_csv(f"{V21}/model_logs/trained_model{suf}.predictions.tsv.gz", sep="\t",
                     usecols=["meta_sample", "taxon", "p_species"])
    tr = tr.merge(pr, on=["meta_sample", "taxon"], how="inner")
    tr = tr[(tr["allele_copy_share"] > 0.2)]
    tr["cls"] = gclass(tr)
    pos = tr[(tr["truth"] == 1) & (tr["cls"] == "real")]
    neg = tr[(tr["truth"] == 0) & (pd.to_numeric(tr["meta_novel_congener"], errors="coerce") == 1)]
    hard_pos = pos[pos["p_species"] < 0.9]
    hard_neg = neg[neg["p_species"] > 0.1]
    print(f"-- {t} v21 training rows on allele copies: real strains {len(pos)} (p<0.9: {len(hard_pos)}), absent beside "
          f"a held-out congener {len(neg)} (p>0.1: {len(hard_neg)}); AUC real strain over FP-like, all / hard:")
    for f in FEATS:
        print(f"   {f:26s} {auc(pos[f], neg[f]):.3f} / {auc(hard_pos[f], hard_neg[f]):.3f}")
