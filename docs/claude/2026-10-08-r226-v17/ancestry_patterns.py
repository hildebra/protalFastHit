#!/usr/bin/env python3
"""The r226 v17 models' errors and the ancestry features: the error budget at the knob by class of row, the new
features (ancestry sites, congener gaps, foreign copies) by class and inside the identity band where strains and
novel congeners overlap, what the false positives and false negatives look like on them, and how v17's rows differ
from v15's (the reads' identity).

    python3 ancestry_patterns.py --build local/v17 [--v15 local/v15] > ancestry_patterns.txt
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
SUFFIX = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
NEW = ["ancestry_sites_per_record", "ancestry_agreement", "ancestry_congener_share", "foreign_scanned_share",
       "foreign_copy_share", "foreign_genus_copy_share", "gap_informative_share", "gap_within_min_share",
       "gap_within_median_share", "gap_position", "untried_candidate_rate", "db_nearest_congener"]
OLD = ["identity", "top_identity", "lu_per_kb", "em_own_share", "excess_scaled_median", "low_mapq_share",
       "cluster_genomes_log10", "third_position_share", "fragments"]
BASE = ["identity", "top_identity", "fragments", "lu_per_kb", "lu_gene_rate3", "excess_scaled_median", "low_mapq_share",
        "em_own_share", "sample_log_fragments", "cluster_genomes_log10", "mean_mapq", "hit_gene_fraction",
        "genus_spill", "relative_distance"]
BINS = [0.95, 0.96, 0.97, 0.98, 0.985]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--build", required=True)
    p.add_argument("--v15", help="the v15 build's folder, for the identity comparison")
    p.add_argument("--types", default="pe,se,pb,ont")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


def auc(y, x):
    y = np.asarray(y).astype(int)
    if len(set(y)) < 2:
        return float("nan")
    return roc_auc_score(y, x)


def load(build, rt):
    """Every row of the training and test tables with its call (model_logs/trained_model<sfx>.calls.tsv.gz: the
    training rows with species held out, the test rows by the final model, at the build's knob)."""
    frames = []
    for which in ("training", "test"):
        path = os.path.join(build, which, TABLES[rt])
        if os.path.isfile(path):
            t = pd.read_csv(path, sep="\t")
            t["set"] = which
            frames.append(t)
    t = pd.concat(frames, ignore_index=True)
    t["taxon"] = t["taxon"].astype(str)
    calls = pd.read_csv(os.path.join(build, "model_logs", f"trained_model{SUFFIX[rt]}.calls.tsv.gz"), sep="\t",
                        usecols=["meta_sample", "taxon", "set", "p", "knob", "call"])
    calls["taxon"] = calls["taxon"].astype(str)
    t = t.merge(calls, on=["set", "meta_sample", "taxon"], how="inner")
    t["truth"] = pd.to_numeric(t["truth"], errors="coerce").fillna(0).astype(int)
    t["call"] = t["call"].astype(int) == 1
    for c in NEW + OLD + BASE:
        if c in t:
            t[c] = pd.to_numeric(t[c], errors="coerce")
    return t


def classes(t):
    present = t["truth"] == 1
    insilico = pd.to_numeric(t["meta_insilico_strain"], errors="coerce").fillna(0) == 1
    rep = pd.to_numeric(t["meta_rep_genome"], errors="coerce").fillna(0) == 1
    level = t["meta_novel_level"].fillna("").astype(str)
    rank = t["meta_relative_rank"].fillna("").astype(str)
    return np.select([present & rep, present & insilico, present, (level == "species") & (rank == "genus"), level != ""],
                     ["present rep", "present insilico", "present strain", "absent novel congener",
                      "absent near novel"], "absent")


def f1(tp, fp, fn):
    return 2 * tp / max(1, 2 * tp + fp + fn)


def budget(t):
    for which in ("training", "test"):
        c = t[t["set"] == which]
        y = c["truth"] == 1
        tp, fp, fn = int((c["call"] & y).sum()), int((c["call"] & ~y).sum()), int((~c["call"] & y).sum())
        print(f"  {which} rows ({'species held out' if which == 'training' else 'final model'}, knob {c['knob'].iloc[0]}): "
              f"{len(c)} rows, {int(y.sum())} present: TP {tp}, FP {fp}, FN {fn}, F1 {f1(tp, fp, fn):.4f}")
        print(f"    FP by class {c.loc[c['call'] & ~y, 'class'].value_counts().to_dict()}; FN by class "
              f"{c.loc[~c['call'] & y, 'class'].value_counts().to_dict()}")
        for name, sel in (("FP", c["call"] & ~y), ("FN", ~c["call"] & y)):
            x = c.loc[sel]
            if len(x):
                print(f"    {name}: fragments median {x['fragments'].median():.0f} (1: {np.mean(x['fragments'] <= 1):.2f}, "
                      f">= 10: {np.mean(x['fragments'] >= 10):.2f}), identity median {x['identity'].median():.4f}, "
                      f"by scenario {x['meta_scenario'].fillna('design').replace('', 'design').value_counts().to_dict()}")


def by_class(t, cols, min_fragments):
    sub = t[t["fragments"] >= min_fragments]
    n = sub["class"].value_counts()
    print(f"  medians by class, rows with >= {min_fragments} fragments ({n.to_dict()}):")
    print(sub.groupby("class")[cols].median().T.round(4).to_string())


def band(t, cols, seed):
    sub = t[t["identity"].between(BINS[0], BINS[-1]) & (t["fragments"] >= 3)
            & t["class"].isin(["present strain", "present insilico", "present rep", "absent novel congener"])].copy()
    print(f"  band {BINS[0]}-{BINS[-1]}, >= 3 fragments: {len(sub)} rows {sub['class'].value_counts().to_dict()}; of the FP "
          f"{int((t['call'] & (t['truth'] == 0) & t['identity'].between(BINS[0], BINS[-1])).sum())} of "
          f"{int((t['call'] & (t['truth'] == 0)).sum())}, of the FN "
          f"{int((~t['call'] & (t['truth'] == 1) & t['identity'].between(BINS[0], BINS[-1])).sum())} of "
          f"{int((~t['call'] & (t['truth'] == 1)).sum())}")
    rows = []
    for col in cols:
        if col not in sub:
            continue
        r = {"feature": col, "whole band (present high)": auc(sub["truth"], sub[col])}
        real = sub[sub["class"] != "present insilico"]
        r["real strains only"] = auc(real["truth"], real[col])
        for lo, hi in zip(BINS[:-1], BINS[1:]):
            b = real[real["identity"].between(lo, hi)]
            r[f"{lo}-{hi}"] = auc(b["truth"], b[col])
        rows.append(r)
    print("  AUC of each feature for present (strains, reps) against absent novel congener, inside the band and in narrow "
          "identity bins (real strains only):")
    print(pd.DataFrame(rows).round(3).to_string(index=False))
    # Does the ancestry group add to the rest, inside the band? Boosted classifier, samples grouped.
    y = sub["truth"].to_numpy()
    g = sub["meta_sample"].to_numpy()
    base = [c for c in BASE if c in sub]
    anc = ["ancestry_sites_per_record", "ancestry_agreement", "ancestry_congener_share"]
    foreign = ["foreign_scanned_share", "foreign_copy_share", "foreign_genus_copy_share"]
    for label, cols_ in (("base", base), ("base + ancestry", base + anc), ("base + foreign", base + foreign),
                         ("base + ancestry + foreign", base + anc + foreign)):
        X = sub[cols_].to_numpy(dtype=float)
        X = np.where(np.isfinite(X), X, np.nan)
        s = np.zeros(len(sub))
        for tr, te in GroupKFold(n_splits=5).split(X, y, g):
            m = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.1, max_leaf_nodes=31, random_state=seed)
            s[te] = m.fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        print(f"    {label}: AUC {roc_auc_score(y, s):.4f}, AP {average_precision_score(y, s):.4f}")


def errors(t, cols):
    t = t.copy()
    y = t["truth"] == 1
    t["error"] = np.select([t["call"] & ~y, ~t["call"] & y, t["call"]], ["FP", "FN", "TP"], "TN")
    sub = t[t["fragments"] >= 3]
    print("  medians of the errors against the right calls, rows with >= 3 fragments:")
    print(sub.groupby("error")[cols].median().T.round(4).to_string())
    e = sub[sub["error"].isin(["FP", "FN"])]
    rows = [{"feature": c, "AUC FN vs FP (FN high)": auc((e["error"] == "FN").astype(int), e[c])} for c in cols if c in e]
    print(pd.DataFrame(rows).round(3).to_string(index=False))


def sites_availability(t):
    sub = t[t["fragments"] >= 3]
    has = sub["ancestry_sites_per_record"] > 0
    print(f"  rows with ancestry sites (>= 3 fragments): {has.mean():.3f}; by class "
          f"{sub.groupby('class')['ancestry_sites_per_record'].apply(lambda x: (x > 0).mean()).round(3).to_dict()}")
    print(f"  sites per record where any: median {sub.loc[has, 'ancestry_sites_per_record'].median():.2f}; agreement "
          f"by class where any: {sub[has].groupby('class')['ancestry_agreement'].median().round(3).to_dict()}; congener share "
          f"{sub[has].groupby('class')['ancestry_congener_share'].median().round(3).to_dict()}")
    d = sub.loc[has, "db_nearest_congener"]
    sub2 = sub[has].copy()
    sub2["near"] = pd.cut(sub2["db_nearest_congener"], [-0.001, 0.02, 0.05, 0.1, 0.15, 1.01])
    print("  agreement by the nearest kept congener's distance and class (median):")
    print(sub2.groupby(["near", "class"], observed=True)["ancestry_agreement"].median().unstack("class").round(3).to_string())


def identity_shift(build, v15, rt):
    rows = []
    for name, b in (("v17", build), ("v15", v15)):
        frames = []
        for which in ("training", "test"):
            path = os.path.join(b, which, TABLES[rt])
            if os.path.isfile(path):
                frames.append(pd.read_csv(path, sep="\t", usecols=["truth", "meta_rep_genome", "meta_insilico_strain",
                                                                   "meta_novel_level", "meta_relative_rank", "identity",
                                                                   "fragments"]))
        x = pd.concat(frames, ignore_index=True)
        x["class"] = classes(x)
        x = x[x["fragments"] >= 3]
        for cls in ("present rep", "present strain", "present insilico", "absent novel congener"):
            s = x.loc[x["class"] == cls, "identity"]
            rows.append({"build": name, "class": cls, "rows": len(s), "10%": s.quantile(0.1), "median": s.median(),
                         "90%": s.quantile(0.9), "0.95-0.985": s.between(0.95, 0.985).mean()})
    print(pd.DataFrame(rows).round(4).to_string(index=False))


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    for rt in opts.types.split(","):
        t = load(opts.build, rt)
        t["class"] = classes(t)
        new = [c for c in NEW if c in t]
        print(f"\n# {rt}: {len(t)} rows, {t['meta_sample'].nunique()} samples, knob {t['knob'].iloc[0]}\n", flush=True)
        print("## The error budget at the knob\n")
        budget(t)
        print("\n## The new features by class\n")
        by_class(t, new + ["identity"], 5)
        print("\n## Inside the identity band 0.95-0.985\n")
        band(t, new + OLD, opts.seed)
        print("\n## The errors on the new features\n")
        errors(t, new + ["identity", "lu_per_kb", "excess_scaled_median"])
        print("\n## Where the ancestry sites exist\n")
        sites_availability(t)
        if opts.v15 and os.path.isfile(os.path.join(opts.v15, "training", TABLES[rt])):
            print("\n## The reads' identity, v17 against v15 (rows with >= 3 fragments)\n")
            identity_shift(opts.build, opts.v15, rt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
