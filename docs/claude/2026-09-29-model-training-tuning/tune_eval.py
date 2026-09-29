#!/usr/bin/env python3
"""Score every design's model (and the shipped one) with protal on the independent test samples and compare.
usage: tune_eval.py TUNE_DIR PROTAL TEST_SET MODEL_NAMES..."""
import glob
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, log_loss

T, PROTAL, TEST, names = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
SHIPPED = "/mnt/c/Users/hildebra/Documents/locDev/protal/scripts/random_forest.xml"
DB = os.path.join(T, "A", "protal_db")


def samples():
    out = []
    for meta in sorted(glob.glob(os.path.join(T, TEST, "points", "*", "sim", "protal.meta"))):
        output_dir = None
        for line in open(meta):
            f = line.rstrip("\n").split("\t")
            if f[0] == "#OUTPUT_DIR":
                output_dir = f[1]
            elif f[0] == "#SAMPLEID":
                header = [x.lstrip("#") for x in f]
            elif not line.startswith("#") and line.strip():
                row = dict(zip(header, f))
                out.append((row["SAMPLEID"], os.path.join(output_dir, "alignments", row["SAM"]), row["PROFILE_TRUTH"]))
    return out


def score(name, model):
    out = os.path.join(T, "score", TEST, name)
    tests = samples()
    if len(glob.glob(os.path.join(out, "*.truth_annotated"))) != len(tests):
        os.makedirs(out, exist_ok=True)
        cmd = [PROTAL, "--db", DB, "--model", model, "--profile_only", ",".join(s[1] for s in tests),
               "--profile_truth", ",".join(s[2] for s in tests), "-o", out, "-t", "8", "--no_strains", "--no_qcmsa"]
        with open(out + ".log", "w") as log:
            if subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode:
                sys.exit(f"protal failed for {name}; see {out}.log")
    parts = []
    for sample, _, _ in tests:
        d = pd.read_csv(os.path.join(out, sample + ".profile.truth_annotated"), sep="\t", float_precision="round_trip")
        d["meta_sample"] = sample
        parts.append(d[["meta_sample", "taxon", "truth", "probability"]])
    return pd.concat(parts, ignore_index=True)


meta = pd.read_csv(os.path.join(T, TEST, "training_data.tsv"), sep="\t", float_precision="round_trip",
                   usecols=lambda c: c.startswith("meta_") or c in ("taxon", "taxon_name", "fragments", "identity"))


def summarize(d, knob):
    y, p = d.truth.to_numpy().astype(int), d.probability.to_numpy()
    call = p >= knob
    tp, fp, fn = int((call & (y == 1)).sum()), int((call & (y == 0)).sum()), int((~call & (y == 1)).sum())
    pres = y == 1
    strain = pres & (d.meta_rep_genome == 0).to_numpy()
    rep = pres & (d.meta_rep_genome == 1).to_numpy()
    congener = (y == 0) & (d.meta_novel_congener == 1).to_numpy()
    other_abs = (y == 0) & (d.meta_novel_congener != 1).to_numpy()
    arch = pres & (d.meta_domain == "Archaea").to_numpy()
    return {"F1": 2 * tp / (2 * tp + fp + fn), "sens": tp / (tp + fn), "prec": tp / max(1, tp + fp),
            "FP/sample": fp / d.meta_sample.nunique(), "sens_strain": call[strain].mean(), "sens_rep": call[rep].mean(),
            "sens_arch": call[arch].mean(), "FP_congener_%": 100 * call[congener].mean(),
            "FP_other_%": 100 * call[other_abs].mean(), "AP": average_precision_score(y, p),
            "logloss": log_loss(y, np.clip(p, 1e-6, 1 - 1e-6), labels=[0, 1])}


def best_knob(d):
    y, p = d.truth.to_numpy().astype(int), d.probability.to_numpy()
    best = (0, 0.5)
    for t in np.round(np.arange(0.05, 0.96, 0.01), 2):
        call = p >= t
        tp, fp, fn = (call & (y == 1)).sum(), (call & (y == 0)).sum(), (~call & (y == 1)).sum()
        best = max(best, (2 * tp / (2 * tp + fp + fn), t))
    return best[1]


rows, depth_rows = [], []
models = [("shipped", SHIPPED)] + [(n, os.path.join(T, n, "trained_model.xml")) for n in names]
for name, model in models:
    if not os.path.exists(model):
        continue
    d = score(name, model).merge(meta, on=["meta_sample", "taxon"], how="left")
    info = {}
    mj = os.path.join(T, name, "trained_model.metrics.json")
    if os.path.exists(mj):
        m = json.load(open(mj))
        info = {"own_F1": m["evaluation"]["species"]["F1"], "own_knob": m["threshold"]["best_F1"]["threshold"],
                "train_taxa": m["data"]["taxa"], "MB": m["model"]["bytes"] / 1e6}
    row = {"model": name, **info, **summarize(d, 0.5)}
    row["test_knob"] = best_knob(d)
    if "own_knob" in info:
        row["F1_at_own_knob"] = summarize(d, info["own_knob"])["F1"]
    row["F1_at_test_knob"] = summarize(d, row["test_knob"])["F1"]
    rows.append(row)
    for pairs, g in d[d.truth == 1].groupby("meta_read_pairs"):
        depth_rows.append({"model": name, "pairs": pairs, "found": f"{int((g.probability >= 0.5).sum())}/{len(g)}"})
    for pairs, g in d[d.truth == 0].groupby("meta_read_pairs"):
        depth_rows.append({"model": name, "pairs": pairs, "FP": int((g.probability >= 0.5).sum())})

ref = meta
print(f"test: {ref.meta_sample.nunique()} samples, {len(ref)} taxa ({int((ref.meta_rep_genome == 0).sum())} present from "
      f"other strains, {int((ref.meta_rep_genome == 1).sum())} from representatives; "
      f"{int((ref.meta_novel_congener == 1).sum())} rows congeneric with an unknown species in the sample)")
pd.set_option("display.width", 250)
print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.4f}", na_rep="-"))
dr = pd.DataFrame(depth_rows)
print("\npresent found at knob 0.5, by read pairs:")
print(dr.dropna(subset=["found"]).pivot(index="pairs", columns="model", values="found").to_string())
print("\nfalse positives at knob 0.5, by read pairs:")
print(dr.dropna(subset=["FP"]).pivot(index="pairs", columns="model", values="FP").to_string())
