#!/usr/bin/env python3
"""Compare the models of eval.sh on the independent test samples, as protal scored them.
usage: eval_models.py EVAL_DIR"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, log_loss

W = sys.argv[1]
WORLDS = {
    "toy": {"genomes": os.path.expanduser("~/audit3/gtdb24/simulation/genomes.tsv"), "train_genera": None},
    "hard": {"genomes": os.path.expanduser("~/audit4/genomes_test.tsv"),
             "train_genera": os.path.expanduser("~/audit4/train_genera.txt")},
}
MODELS = ["shipped", "old_gtdb", "old_default", "new", "new_all"]
LABEL = {"shipped": "shipped model (2024)", "old_gtdb": "old trainer, GTDB workflow (normalized, 512 trees)",
         "old_default": "old trainer, defaults (all features, 256 trees)", "new": "new trainer (normalized, 64 trees)",
         "new_all": "new trainer, --features all"}


def domains(path):
    """species name -> domain, from the lineages of a genome table"""
    out = {}
    for line in open(path):
        lineage = next((f for f in line.rstrip("\n").split("\t") if f.startswith("d__") and ";s__" in f), None)
        if lineage:
            out[lineage.split(";")[-1]] = lineage.split(";")[0][3:]
    return out


def self_report(world, model):
    if model.startswith("old"):
        log = open(f"{W}/{world}/{model}.log").read()
        m = re.search(r"Test\s+sens=([\d.]+) prec=([\d.]+) f1=([\d.]+)", log)
        return float(m.group(3)) if m else None
    if model.startswith("new"):
        path = f"{W}/{world}/{model}.metrics.json"
        if os.path.exists(path):
            return json.load(open(path))["evaluation"]["species"]["F1"]
    return None


def load(world, model):
    files = glob.glob(f"{W}/{world}/score_{model}/**/*.truth_annotated", recursive=True)
    if not files:
        return None
    parts = []
    for f in files:
        d = pd.read_csv(f, sep="\t", float_precision="round_trip")
        point = f.split(f"score_{model}/")[1].split("/")[0]
        d["point"] = point
        d["pairs"] = int(point.split("_p")[1])
        d["sample"] = os.path.basename(f).split(".profile")[0]
        parts.append(d)
    return pd.concat(parts, ignore_index=True)


def summarize(d, knob=0.5):
    y, p = d.truth.to_numpy().astype(int), d.probability.to_numpy()
    call = p >= knob
    tp, fp, fn = int((call & (y == 1)).sum()), int((call & (y == 0)).sum()), int((~call & (y == 1)).sum())
    arch = (d.domain == "Archaea").to_numpy() & (y == 1)
    bac = (d.domain == "Bacteria").to_numpy() & (y == 1)
    return {"F1": 2 * tp / (2 * tp + fp + fn), "sens": tp / (tp + fn), "prec": tp / (tp + fp) if tp + fp else np.nan,
            "FP": fp, "FP/sample": fp / d["sample"].nunique(), "sens_arch": call[arch].mean(), "sens_bac": call[bac].mean(),
            "AP": average_precision_score(y, p), "logloss": log_loss(y, np.clip(p, 1e-6, 1 - 1e-6), labels=[0, 1])}


for world, cfg in WORLDS.items():
    dom = domains(cfg["genomes"])
    genera = set(open(cfg["train_genera"]).read().split()) if cfg["train_genera"] else None
    runtime = {}
    if os.path.exists(f"{W}/{world}/runtime.txt"):
        for line in open(f"{W}/{world}/runtime.txt"):
            model, size, secs, kb = line.split()
            runtime[model] = (int(size) / 1e6, float(secs), int(kb) / 1e3)
    rows, by_seen, by_depth = [], [], []
    data = {}
    for model in MODELS:
        d = load(world, model)
        if d is None:
            continue
        d["domain"] = d.taxon_name.map(dom).fillna("unknown")
        data[model] = d
        timef = f"{W}/{world}/{model}.time"
        row = {"model": LABEL[model]}
        if os.path.exists(timef):
            row["train_s"] = float(open(timef).read().split()[-1])
        if os.path.exists(f"{W}/{world}/new_fit.time") and model == "new":
            row["fit_only_s"] = float(open(f"{W}/{world}/new_fit.time").read().split()[-1])
        if model in runtime:
            row["MB"], row["protal_s"], row["protal_MB_RAM"] = runtime[model]
        row["self_F1"] = self_report(world, model)
        row.update(summarize(d))
        rows.append(row)
        if genera is not None:
            d["genus_seen"] = d.taxon_name.str.slice(3).str.split(" ").str[0].isin(genera)
            for seen, g in d.groupby("genus_seen"):
                s = summarize(g)
                by_seen.append({"model": model, "genus in training": seen, "present": int(g.truth.sum()),
                                "absent": int((g.truth == 0).sum()), "sens": s["sens"], "FP": s["FP"],
                                "sens_arch": s["sens_arch"], "sens_bac": s["sens_bac"]})
        for pairs, g in d[d.truth == 1].groupby("pairs"):
            by_depth.append({"model": model, "pairs": pairs,
                             "archaea found": f"{int((g[g.domain == 'Archaea'].probability >= 0.5).sum())}/{int((g.domain == 'Archaea').sum())}",
                             "bacteria found": f"{int((g[g.domain == 'Bacteria'].probability >= 0.5).sum())}/{int((g.domain == 'Bacteria').sum())}"})
    if not data:
        continue
    ref = next(iter(data.values()))
    print(f"\n### {world} world: {ref['sample'].nunique()} test samples, {len(ref)} taxa "
          f"({int(ref.truth.sum())} present, {int((ref.truth == 0).sum())} absent)")
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.4f}", na_rep="-"))
    if by_seen:
        print("\nby whether the genus was in the training samples:")
        print(pd.DataFrame(by_seen).to_string(index=False, float_format=lambda v: f"{v:.4f}", na_rep="-"))
    print("\npresent taxa found, by depth:")
    t = pd.DataFrame(by_depth)
    print(t.pivot(index="pairs", columns="model", values="archaea found").to_string())
    print(t.pivot(index="pairs", columns="model", values="bacteria found").to_string())
    # the same taxa in every model's scoring?
    keys = {m: set(zip(d["sample"], d.taxon)) for m, d in data.items()}
    same = all(k == keys[next(iter(keys))] for k in keys.values())
    print(f"\nsame taxa scored by every model: {same}")
