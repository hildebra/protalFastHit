"""v23's errors with species held out (the build's own models, model_logs/trained_model*.calls.tsv.gz, knob 0.5): by
kind of sample (design or scenario), by the taxon's fragments, and by the distance to the nearest held-out congener.

usage: error_anatomy.py <build dir>
"""
import sys
import numpy as np
import pandas as pd

B = sys.argv[1]
pd.set_option("display.width", 250)
for t in ("pe", "se", "pb", "ont"):
    suf = "" if t == "pe" else "_" + t
    c = pd.read_csv(f"{B}/model_logs/trained_model{suf}.calls.tsv.gz", sep="\t",
                    usecols=["meta_sample", "meta_scenario", "meta_rep_genome", "meta_insilico_strain",
                             "meta_novel_congener", "meta_novel_distance", "taxon_name", "truth", "set", "p"],
                    dtype={"meta_scenario": str})
    c = c[c["set"] == "training"]
    frag = pd.read_csv(f"{B}/work/training/training_data{suf}.tsv", sep="\t",
                       usecols=["meta_sample", "taxon_name", "fragments"])
    c = c.merge(frag.drop_duplicates(["meta_sample", "taxon_name"]), on=["meta_sample", "taxon_name"], how="left")
    call = c["p"] >= 0.5
    present = c["truth"] == 1
    c["kind"] = np.where(c["meta_scenario"].fillna("") == "", "design", "scenario:" + c["meta_scenario"].fillna(""))
    c["err"] = np.where(present, ~call, call)
    c["type"] = np.where(present, "FN", "FP")
    c["frag"] = pd.cut(c["fragments"], [-1, 1, 3, 9, 99, 1e12], labels=["1", "2-3", "4-9", "10-99", ">=100"])
    print(f"== {t}: errors with species held out by kind of sample (samples, FP, FN, taxa)")
    k = c.groupby("kind").agg(samples=("meta_sample", "nunique"), taxa=("err", "size"),
                              FP=("err", lambda e: int((e & (c.loc[e.index, "type"] == "FP")).sum())),
                              FN=("err", lambda e: int((e & (c.loc[e.index, "type"] == "FN")).sum())))
    k["FP share"] = (k["FP"] / k["FP"].sum()).round(3)
    k["FN share"] = (k["FN"] / k["FN"].sum()).round(3)
    print(k.to_string())
    print(f"-- {t}: by the taxon's fragments: FP among absent, FN among present (count, rate)")
    rows = []
    for f, g in c.groupby("frag", observed=True):
        a, p = g[g["type"] == "FP"], g[g["type"] == "FN"]
        rows.append({"fragments": f, "absent": len(a), "FP": int(a["err"].sum()), "FP rate": round(a["err"].mean(), 4),
                     "present": len(p), "FN": int(p["err"].sum()), "FN rate": round(p["err"].mean(), 4)})
    print(pd.DataFrame(rows).to_string(index=False))
    a = c[(c["type"] == "FP")]
    d = pd.to_numeric(a["meta_novel_distance"], errors="coerce")
    a = a.assign(dist=pd.cut(d, [0, 0.02, 0.05, 0.15, 10], labels=["<0.02", "0.02-0.05", "0.05-0.15", ">0.15"]))
    a["dist"] = a["dist"].cat.add_categories("none").fillna("none")
    print(f"-- {t}: FP by the distance to the nearest held-out congener in the sample")
    print(a.groupby("dist", observed=True)["err"].agg(["size", "sum", "mean"]).rename(
        columns={"size": "absent", "sum": "FP", "mean": "rate"}).round(4).to_string())
    print()
