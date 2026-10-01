#!/usr/bin/env python3
"""The new build's features, computed by protal itself (profile_new_build.sh, tables_new_build.py): F1, precision and
sensitivity of the forest with and without them, cross-validated on the training table (species held out) and on the
independent test set, over forest seeds; the V2 build (old binary, shipped features) for reference.

    python3 run_new_build.py NEW_OUT [seeds]
"""
import os
import sys

import numpy as np
import pandas as pd

import exp_lib as L

NEW = os.path.expanduser(sys.argv[1])
SEEDS = int(sys.argv[2]) if len(sys.argv) > 2 else 5
BASE = list(L.NORMALIZED_FEATURES)
VARIANTS = {
    "shipped features": [],
    "+ MAPQ": ["mean_mapq", "low_mapq_share"],
    "+ congener fit": ["congener_fit_share", "other_genus_fit_share"],
    "+ MAPQ, congener fit": ["mean_mapq", "low_mapq_share", "congener_fit_share", "other_genus_fit_share"],
    "+ linked reads": ["linked_share"],
    "+ all five": ["mean_mapq", "low_mapq_share", "congener_fit_share", "other_genus_fit_share", "linked_share"],
}
pd.set_option("display.width", 230)


def load(root, split, rt):
    df = pd.read_csv(f"{root}/{split}/{L.TABLE[rt]}", sep="\t", float_precision="round_trip", low_memory=False)
    df["truth"] = df["truth"].astype(str).str.lower().isin(["1", "true"]).astype(int)
    return df


def rates(s):
    return {"precision": s["TP"] / max(1, s["TP"] + s["FP"]), "sensitivity": s["TP"] / max(1, s["TP"] + s["FN"])}


rows = []
for rt in L.READ_TYPES:
    builds = {"new build": (load(NEW, "training", rt), load(NEW, "test", rt))}
    if os.path.exists(f"{L.V2}/training/{L.TABLE[rt]}"):  # V2's own tables, while they exist
        builds = {"V2 build": (load(L.V2, "training", rt), load(L.V2, "test", rt)), **builds}
    for build, (train, test) in builds.items():
        y_tr, y_te = train["truth"].to_numpy(), test["truth"].to_numpy()
        base_calls = {}
        for name, add in VARIANTS.items():
            if build == "V2 build" and add:
                continue
            feats = BASE + add
            for seed in range(1, SEEDS + 1):
                cv = L.scores(y_tr, L.cv_predict(train, feats, seed) >= L.KNOB)
                call = L.test_predict(train, test, feats, seed) >= L.KNOB
                te = L.scores(y_te, call)
                if name == "shipped features":
                    base_calls[seed] = call
                d = L.bootstrap_diff(test, base_calls[seed], call, n=500, seed=seed) if add else (0.0, 0.0, 0.0)
                rows.append({"read_type": rt, "build": build, "features": name, "seed": seed, "cv_F1": cv["F1"],
                             "cv_FP": cv["FP"], "cv_FN": cv["FN"], "test_F1": te["F1"], "test_FP": te["FP"],
                             "test_FN": te["FN"], **{f"test_{k}": v for k, v in rates(te).items()},
                             "dF1": d[0], "lo": d[1], "hi": d[2]})
            r = pd.DataFrame([x for x in rows if x["read_type"] == rt and x["build"] == build and x["features"] == name])
            print(f"{rt} {build:9s} {name:22s} CV F1 {r.cv_F1.mean():.4f}  test F1 {r.test_F1.mean():.4f} "
                  f"precision {r.test_precision.mean():.4f} sensitivity {r.test_sensitivity.mean():.4f} "
                  f"(FP {r.test_FP.mean():.1f}, FN {r.test_FN.mean():.1f})", flush=True)

out = pd.DataFrame(rows)
out.to_csv(f"{NEW}/new_build_features.tsv", sep="\t", index=False)
summary = out.groupby(["read_type", "build", "features"], sort=False).agg(
    cv_F1=("cv_F1", "mean"), cv_FP=("cv_FP", "mean"), cv_FN=("cv_FN", "mean"), test_F1=("test_F1", "mean"),
    test_F1_sd=("test_F1", "std"), precision=("test_precision", "mean"), sensitivity=("test_sensitivity", "mean"),
    FP=("test_FP", "mean"), FN=("test_FN", "mean"), dF1=("dF1", "mean"), lo=("lo", "mean"), hi=("hi", "mean"))
print()
print(summary.round(4).to_string())
