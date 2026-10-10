"""How well each new feature alone separates present real strains from absent species beside a held-out congener (the
classes that carry most FN and FP), on v22's training rows with species held out: AUC over all such rows and over the
hard ones (strains with p < 0.9, absent with p > 0.1 by the build's model, species held out). Also how many rows have
the feature known (not -1). An AUC below 0.5 means the strains score lower.

usage: feature_auc.py <v22 dir> <read type>
"""
import sys
import numpy as np
import pandas as pd

V22, T = sys.argv[1], sys.argv[2]
SUF = "" if T == "pe" else "_" + T
FEATS = ["ancestry_agreement", "ancestry_congener_share", "ancestry_sites_per_record", "ancestry_fixed_gain",
         "ancestry_fixed_agreement", "conserved_mismatch_ratio", "conserved_mismatch_rate", "ancestry_agreement_weighted",
         "nonsynonymous_share", "nonsynonymous_conserved_rate", "allele_explained_share", "column_weight_coverage",
         "identity"]


def auc(pos, neg):
    """P(pos > neg) + 0.5 P(tie)."""
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    ranks = pd.Series(np.concatenate([pos, neg])).rank().to_numpy()
    return (ranks[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


tab = pd.read_csv(f"{V22}/work/training/training_data{SUF}.tsv", sep="\t",
                  usecols=["meta_sample", "taxon", "truth", "meta_rep_genome", "meta_insilico_strain",
                           "meta_novel_congener"] + FEATS)
calls = pd.read_csv(f"{V22}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t",
                    usecols=["meta_sample", "taxon", "set", "p"])
calls = calls[calls["set"] == "training"]
tab = tab.merge(calls[["meta_sample", "taxon", "p"]], on=["meta_sample", "taxon"], how="inner")
rep = pd.to_numeric(tab["meta_rep_genome"], errors="coerce") == 1
ins = pd.to_numeric(tab["meta_insilico_strain"], errors="coerce") == 1
cong = pd.to_numeric(tab["meta_novel_congener"], errors="coerce") == 1
strain = (tab["truth"] == 1) & ~rep & ~ins
beside = (tab["truth"] == 0) & cong
hard_s, hard_b = strain & (tab["p"] < 0.9), beside & (tab["p"] > 0.1)
print(f"== {T}: real strains {strain.sum()} (hard {hard_s.sum()}), absent beside a held-out congener {beside.sum()} "
      f"(hard {hard_b.sum()}); AUC of strains over absent")
rows = []
for f in FEATS:
    v = tab[f].to_numpy(float)
    known = v != -1
    r = {"feature": f, "known strains": f"{known[strain].mean():.2f}", "known absent": f"{known[beside].mean():.2f}",
         "AUC all": round(auc(v[strain & known], v[beside & known]), 3),
         "AUC hard": round(auc(v[hard_s & known], v[hard_b & known]), 3)}
    rows.append(r)
print(pd.DataFrame(rows).to_string(index=False))
