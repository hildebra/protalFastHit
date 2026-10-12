"""How far the simulated strains are from the representatives the database holds: the read identity protal measured
(identity, top_identity) of the present taxa with 11 or more fragments by the genome simulated, and of the absent taxa
beside a held-out congener (meta_novel_congener) that the model called and did not call, per build (training tables).

usage: strain_identity.py <local dir> <read type> <build> [<build> ...]
"""
import glob
import sys
import numpy as np
import pandas as pd

L, T, builds = sys.argv[1], sys.argv[2], sys.argv[3:]
SUF = "" if T == "pe" else "_" + T
Q = (0.05, 0.25, 0.5)
for v in builds:
    d = glob.glob(f"{L}/{v}/protal*_r226_v*")[0]
    t = pd.read_csv(f"{d}/work/training/training_data{SUF}.tsv", sep="\t",
                    usecols=["meta_sample", "taxon_name", "truth", "meta_rep_genome", "meta_insilico_strain",
                             "meta_novel_congener", "fragments", "identity", "top_identity"])
    c = pd.read_csv(f"{d}/model_logs/trained_model{SUF}.calls.tsv.gz", sep="\t",
                    usecols=["meta_sample", "taxon_name", "set", "p"])
    c = c[c["set"] == "training"].drop_duplicates(["meta_sample", "taxon_name"])
    t = t.merge(c[["meta_sample", "taxon_name", "p"]], on=["meta_sample", "taxon_name"], how="left")
    rep = pd.to_numeric(t["meta_rep_genome"], errors="coerce") == 1
    ins = pd.to_numeric(t["meta_insilico_strain"], errors="coerce") == 1
    cong = pd.to_numeric(t["meta_novel_congener"], errors="coerce") == 1
    present = t["truth"] == 1
    deep = t["fragments"] > 10
    groups = {
        "representative, >10 fragments": present & rep & deep,
        "real strain, >10 fragments": present & ~rep & ~ins & deep,
        "in-silico strain, >10 fragments": present & ins & deep,
        "absent beside held-out congener, called": ~present & cong & (t["p"] >= 0.5),
        "absent beside held-out congener, not called, >2 fragments": ~present & cong & (t["p"] < 0.5)
                                                                     & (t["fragments"] > 2),
    }
    rows = []
    for name, m in groups.items():
        g = t[m]
        r = {"class": name, "taxa": len(g)}
        for col in ("identity", "top_identity"):
            for q in Q:
                r[f"{col} {int(q * 100)}%"] = g[col].quantile(q)
        r["identity < 0.98"] = (g["identity"] < 0.98).mean()
        rows.append(r)
    print(f"== {T} {v} (training samples, species held out for the calls)")
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.4f}"))
