#!/usr/bin/env python3
"""New model features against the shipped feature set, per read type: species-held-out cross-validation on the
training table and the independent test set, over forest seeds; paired bootstrap of the test F1 difference.

    python3 run_features.py [seeds] [--with-secondary]
"""
import sys

import numpy as np
import pandas as pd

import exp_lib as L

SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 5
SECONDARY = "--with-secondary" in sys.argv
BASE = list(L.NORMALIZED_FEATURES)
pd.set_option("display.width", 220)


VARIANTS = {
    "base": [],
    "identity_z": ["identity_z"],
    "mapq": ["mean_mapq", "low_mapq_share"],
    "mate_concordance": ["mate_concordance"],
    "clipping": ["clipped_share"],
}
if SECONDARY:
    VARIANTS["congener_fit"] = ["congener_share_d1", "other_genus_share_d1"]
VARIANTS["all"] = sorted({c for v in VARIANTS.values() for c in v})

rows = []
for rt in L.READ_TYPES:
    train, test = L.load("training", rt), L.load("test", rt)
    model = L.identity_model(train)
    for df, split in ((train, "training"), (test, "test")):
        df["identity_z"] = L.identity_z(df, model)
    prim_tr, prim_te = L.primary_features("training", rt), L.primary_features("test", rt)
    train = train.merge(prim_tr, on=["meta_sample", "taxon"], how="left")
    test = test.merge(prim_te, on=["meta_sample", "taxon"], how="left")
    if SECONDARY:
        train = train.merge(L.secondary_features("training", rt), on=["meta_sample", "taxon"], how="left")
        test = test.merge(L.secondary_features("test", rt), on=["meta_sample", "taxon"], how="left")
    extra = sorted({c for v in VARIANTS.values() for c in v})
    for df in (train, test):
        missing = df[extra].isna().any(axis=1).sum()
        if missing:
            print(f"{rt}: {missing} rows without SAM features (filled with 0)")
        df[extra] = df[extra].fillna(0.0)
    print(f"{rt}: identity of present taxa: mu {model[0]:.4f}, sd between taxa {model[1]:.4f}, sd per fragment {model[2]:.4f}; "
          f"training {len(train)} rows, test {len(test)} rows", flush=True)
    y_tr, y_te = train["truth"].to_numpy(), test["truth"].to_numpy()
    base_calls = {}
    for name, add in VARIANTS.items():
        if rt != "pe" and add == ["mate_concordance"]:
            continue
        feats = BASE + [c for c in add if c not in BASE]
        for seed in range(1, SEEDS + 1):
            cv = L.scores(y_tr, L.cv_predict(train, feats, seed) >= L.KNOB)
            p = L.test_predict(train, test, feats, seed)
            call = p >= L.KNOB
            te = L.scores(y_te, call)
            if name == "base":
                base_calls[seed] = call
            diff = L.bootstrap_diff(test, base_calls[seed], call, n=500, seed=seed) if name != "base" else (0, 0, 0)
            rows.append({"read_type": rt, "variant": name, "seed": seed, "cv_F1": cv["F1"], "cv_FP": cv["FP"],
                         "cv_FN": cv["FN"], "test_F1": te["F1"], "test_FP": te["FP"], "test_FN": te["FN"],
                         "test_dF1": diff[0], "test_dF1_lo": diff[1], "test_dF1_hi": diff[2]})
        r = pd.DataFrame([x for x in rows if x["read_type"] == rt and x["variant"] == name])
        print(f"  {name:16s} CV F1 {r.cv_F1.mean():.4f} (FP {r.cv_FP.mean():.0f}, FN {r.cv_FN.mean():.0f})  "
              f"test F1 {r.test_F1.mean():.4f} (FP {r.test_FP.mean():.1f}, FN {r.test_FN.mean():.1f})", flush=True)

out = pd.DataFrame(rows)
out.to_csv(f"{L.EXP}/features_{'with' if SECONDARY else 'without'}_secondary.tsv", sep="\t", index=False)
summary = out.groupby(["read_type", "variant"], sort=False).agg(
    cv_F1=("cv_F1", "mean"), cv_F1_sd=("cv_F1", "std"), cv_FP=("cv_FP", "mean"), cv_FN=("cv_FN", "mean"),
    test_F1=("test_F1", "mean"), test_F1_sd=("test_F1", "std"), test_FP=("test_FP", "mean"), test_FN=("test_FN", "mean"),
    test_dF1=("test_dF1", "mean"), dF1_lo=("test_dF1_lo", "mean"), dF1_hi=("test_dF1_hi", "mean"))
base = summary.xs("base", level="variant")
summary["cv_dF1"] = summary["cv_F1"] - summary.index.get_level_values("read_type").map(base["cv_F1"]).to_numpy()
print()
print(summary.round(4).to_string())
