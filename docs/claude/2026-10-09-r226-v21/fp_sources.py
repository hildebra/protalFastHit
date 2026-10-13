"""Where the FP's reads come from, v19 against v21 (error_reads taxa.tsv.gz): the source species with most fragments
on the FP taxon, as a held-out species (not in the database) or a present species of the database by the genome it
was simulated from (representative, real strain, in-silico strain, from the calls' truth rows of the same sample).
FP per 1000 present taxa of each class gives the spill per class, to separate the strain mix from the rest.

usage: fp_sources.py <local dir> <read type>
"""
import sys
import numpy as np
import pandas as pd

L, T = sys.argv[1], sys.argv[2]
SUF = "" if T == "pe" else "_" + T


def gclass(df):
    rep = pd.to_numeric(df["meta_rep_genome"], errors="coerce") == 1
    ins = pd.to_numeric(df["meta_insilico_strain"], errors="coerce") == 1
    return np.where(rep, "rep", np.where(ins, "insilico", "real"))


for v in ("v19", "v21"):
    base = f"{L}/{v}/protal0.7.9_r226_{v}/model_logs"
    calls = pd.read_csv(f"{base}/trained_model{SUF}.calls.tsv.gz", sep="\t",
                        usecols=["meta_sample", "meta_rep_genome", "meta_insilico_strain", "taxon_name", "truth", "set"])
    pres = calls[calls["truth"] == 1].copy()
    pres["cls"] = gclass(pres)
    cls = dict(zip(zip(pres["meta_sample"], pres["taxon_name"]), pres["cls"]))
    n_cls = pres[pres["set"] == "training"]["cls"].value_counts()
    tx = pd.read_csv(f"{base}/error_reads/{T}/taxa.tsv.gz", sep="\t", dtype=str)
    fp = tx[(tx["error"] == "FP") & (tx["set"] == "training")]
    src = []
    for s, f in zip(fp["sample"], fp["fragments_on_taxon_from"].fillna("")):
        best, bn = None, -1
        for part in f.split(","):
            if ":" not in part:
                continue
            name, n = part.rsplit(":", 1)
            if int(n) > bn:
                best, bn = name, int(n)
        if best is None:
            src.append("none")
        elif best.endswith("(not in the database)"):
            src.append("held-out species")
        else:
            src.append("present " + cls.get((s, best), "absent?"))
    c = pd.Series(src).value_counts()
    print(f"== {T} {v}: training FP (species held out, build knob) {len(fp)} by the source of most of their reads:")
    for k, n in c.items():
        per = ""
        if k.startswith("present ") and k[8:] in n_cls:
            per = f", {1000 * n / n_cls[k[8:]]:.1f} per 1000 present {k[8:]}"
        print(f"   {k:22s} {n:6d}{per}")
    print("   present taxa: " + ", ".join(f"{k} {n}" for k, n in n_cls.items()))
