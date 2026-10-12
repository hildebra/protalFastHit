"""The species' depth error (depth / true_depth - 1, model_logs/composition_species*.tsv.gz) of the present, called
species by the genome simulated (representative, real strain, in-silico strain; meta_* of the calls), per build.

usage: depth_by_class.py <local dir> <read type> <build> [<build> ...]
"""
import glob
import sys
import pandas as pd

L, T, builds = sys.argv[1], sys.argv[2], sys.argv[3:]
SUF = "" if T == "pe" else "_" + T
for v in builds:
    d = glob.glob(f"{L}/{v}/protal*_r226_v*")[0]
    comp = pd.read_csv(f"{d}/model_logs/composition_species{SUF}.tsv.gz", sep="\t",
                       usecols=["set", "scenario", "sample", "species", "present", "called", "true_depth", "depth"],
                       dtype={"scenario": str})
    comp = comp[(comp["present"] == 1) & (comp["called"] == 1) & comp["depth"].notna() & (comp["true_depth"] > 0)]
    calls = pd.read_csv(f"{d}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t",
                        usecols=["meta_sample", "taxon_name", "meta_rep_genome", "meta_insilico_strain", "truth"])
    calls = calls[calls["truth"] == 1].drop_duplicates(["meta_sample", "taxon_name"])
    m = comp.merge(calls, left_on=["sample", "species"], right_on=["meta_sample", "taxon_name"], how="inner")
    rep = pd.to_numeric(m["meta_rep_genome"], errors="coerce") == 1
    ins = pd.to_numeric(m["meta_insilico_strain"], errors="coerce") == 1
    m["class"] = "real strain"
    m.loc[rep, "class"] = "representative"
    m.loc[ins, "class"] = "in-silico strain"
    m["err"] = m["depth"] / m["true_depth"] - 1
    m = m[m["scenario"].fillna("") == ""]
    print(f"== {T} {v}: species' depth / true depth - 1 of the present, called species (design samples), "
          f"{len(m)} of {len(comp[comp['scenario'].fillna('') == ''])} matched")
    rows = []
    for (s, c), g in m.groupby(["set", "class"]):
        rows.append({"set": s, "class": c, "species": len(g), "share": len(g) / (m["set"] == s).sum(),
                     "median": g["err"].median(), "10%": g["err"].quantile(0.1), "90%": g["err"].quantile(0.9)})
    for s, g in m.groupby("set"):
        rows.append({"set": s, "class": "all", "species": len(g), "share": 1.0, "median": g["err"].median(),
                     "10%": g["err"].quantile(0.1), "90%": g["err"].quantile(0.9)})
    print(pd.DataFrame(rows).sort_values(["set", "class"]).to_string(index=False, float_format=lambda x: f"{x:+.4f}"))
