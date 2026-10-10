"""v21 against v22, the builds' own models (model_logs/trained_model*.calls.tsv.gz) at knob 0.5: the present taxa by the
genome simulated and the absent ones by their neighbourhood (as in variants.py), their counts and error rates, with
species held out and on the test set (design samples).

usage: builds_by_class.py <local dir> <read type>
"""
import sys
import numpy as np
import pandas as pd

L, T = sys.argv[1], sys.argv[2]
SUF = "" if T == "pe" else "_" + T
KNOB = 0.5


def classes(c, has):
    present = c["truth"] == 1
    rep = pd.to_numeric(c["meta_rep_genome"], errors="coerce") == 1
    ins = pd.to_numeric(c["meta_insilico_strain"], errors="coerce") == 1
    cong = pd.to_numeric(c["meta_novel_congener"], errors="coerce") == 1
    genus_rel = c["meta_relative_rank"].fillna("") == "genus"
    alle = c["taxon_name"].map(has).fillna(0) > 0.2
    return {
        "FN real strain, alleles": present & ~rep & ~ins & alle,
        "FN real strain, none": present & ~rep & ~ins & ~alle,
        "FN in-silico strain": present & ins,
        "FN representative": present & rep,
        "FP beside held-out congener": ~present & cong,
        "FP beside present congeners only": ~present & ~cong & genus_rel,
        "FP no congener in sample": ~present & ~cong & ~genus_rel,
    }


res = {}
for v in ("v21", "v22"):
    d = f"{L}/{v}/protal0.7.9_r226_{v}"
    tab = pd.concat([pd.read_csv(f"{d}/work/{s}/training_data{SUF}.tsv", sep="\t",
                                 usecols=["taxon_name", "fragments", "allele_copy_share"]) for s in ("training", "test")])
    has = tab[tab["fragments"] >= 3].groupby("taxon_name")["allele_copy_share"].max()
    c = pd.read_csv(f"{d}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t",
                    usecols=["meta_sample", "meta_scenario", "meta_rep_genome", "meta_insilico_strain",
                             "meta_novel_congener", "meta_relative_rank", "taxon_name", "truth", "set", "p"],
                    dtype={"meta_scenario": str, "meta_relative_rank": str})
    sc = c["meta_scenario"].fillna("")
    groups = {"held out": c["set"] == "training", "test": (c["set"] == "test") & (sc == "")}
    call = c["p"].to_numpy() >= KNOB
    for g, gm in groups.items():
        for name, m in classes(c, has).items():
            m = (m & gm).to_numpy()
            err = (~call[m]).sum() if name.startswith("FN") else call[m].sum()
            res[(g, name, v)] = (int(m.sum()), int(err))

for g in ("held out", "test"):
    print(f"== {T}, {g}: class size (share of present or absent) and error rate (count), v21 -> v22")
    rows = []
    for name in [k[1] for k in res if k[0] == g and k[2] == "v21"]:
        r = {"class": name}
        for v in ("v21", "v22"):
            tot = sum(res[(g, k, v)][0] for k in [x[1] for x in res if x[0] == g and x[2] == v]
                      if k[:2] == name[:2])
            n, e = res[(g, name, v)]
            r[f"{v} n"] = f"{n} ({n / max(tot, 1):.1%})"
            r[f"{v} rate"] = f"{e / max(n, 1):.4f} ({e})"
        rows.append(r)
    print(pd.DataFrame(rows).to_string(index=False))

# The mix: each build's error rates applied to the other's class shares (present and absent apart), at the build's
# totals, so that the F1 difference splits into what the rates changed and what the mix did.
print(f"== {T}, held out: F1 at knob 0.5 from the class rates and shares (rates of the row's build, shares of the column's)")
g = "held out"
names = [k[1] for k in res if k[0] == g and k[2] == "v21"]


def f1_from(rate_v, mix_v):
    tot = {p: sum(res[(g, k, mix_v)][0] for k in names if k[:2] == p) for p in ("FN", "FP")}
    fn = sum(tot["FN"] * res[(g, k, mix_v)][0] / tot["FN"] * res[(g, k, rate_v)][1] / max(res[(g, k, rate_v)][0], 1)
             for k in names if k.startswith("FN"))
    fp = sum(tot["FP"] * res[(g, k, mix_v)][0] / tot["FP"] * res[(g, k, rate_v)][1] / max(res[(g, k, rate_v)][0], 1)
             for k in names if k.startswith("FP"))
    tp = tot["FN"] - fn
    return 2 * tp / (2 * tp + fp + fn)


for rv in ("v21", "v22"):
    print(f"   rates {rv}: " + "  ".join(f"mix {mv} {f1_from(rv, mv):.4f}" for mv in ("v21", "v22")))
