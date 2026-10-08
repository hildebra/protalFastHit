#!/usr/bin/env python3
"""The r226 v18 models and the ancestry features (consensus sites, comparison chained across indels, indel sites; v17
had the nearest congener's sites on the main diagonal only): the features by class of row in v18 and v17, their AUC
inside the identity band where strains and novel congeners overlap, what a boosted classifier inside the band gains
from them (v17's three, v18's five, the leak-free foreign rates), the errors on them, and the errors by scenario. Uses
the v17 report's ancestry_patterns.py for loading the tables with the models' calls and classing the rows.

    python3 ancestry_v18.py --build local/v18 --v17 local/v17 > ancestry_v18.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "2026-10-08-r226-v17"))
import ancestry_patterns as ap  # noqa: E402

ANC3 = ["ancestry_sites_per_record", "ancestry_agreement", "ancestry_congener_share"]
INDEL = ["ancestry_indel_sites_per_record", "ancestry_indel_congener_share"]
ANC5 = ANC3 + INDEL
FOREIGN = ["foreign_scanned_share", "foreign_copy_share", "foreign_genus_copy_share"]
GAPS = ["gap_position", "gap_within_min_share", "untried_candidate_rate"]
OTHER = ["identity", "top_identity", "lu_per_kb", "em_own_share", "excess_scaled_median", "db_nearest_congener",
         "db_congeners_01", "fragments"]
CLASSES = ["present rep", "present strain", "present insilico", "absent novel congener", "absent near novel", "absent"]


def numeric(t, cols):
    for c in cols:
        if c in t:
            t[c] = pd.to_numeric(t[c], errors="coerce")
    return t


def load(build, rt):
    t = ap.load(build, rt)
    t["class"] = ap.classes(t)
    return numeric(t, ANC5 + FOREIGN + GAPS + OTHER)


def by_class(t, name, cols, min_fragments=5):
    sub = t[t["fragments"] >= min_fragments]
    cols = [c for c in cols if c in sub]
    m = sub.groupby("class")[cols].median().T.reindex(columns=[c for c in CLASSES if c in set(sub["class"])])
    print(f"  {name}: medians by class, rows with >= {min_fragments} fragments "
          f"({sub['class'].value_counts().reindex(CLASSES).dropna().astype(int).to_dict()})")
    print(m.round(4).to_string())


def band_rows(t):
    return t[t["identity"].between(ap.BINS[0], ap.BINS[-1]) & (t["fragments"] >= 3)
             & t["class"].isin(["present strain", "present insilico", "present rep", "absent novel congener"])].copy()


def band_auc(t, name, cols):
    sub = band_rows(t)
    real = sub[sub["class"] != "present insilico"]
    rows = []
    for col in cols:
        if col not in sub:
            continue
        r = {"feature": col, "build": name, "whole band": ap.auc(sub["truth"], sub[col].fillna(-1)),
             "real strains": ap.auc(real["truth"], real[col].fillna(-1))}
        for lo, hi in zip(ap.BINS[:-1], ap.BINS[1:]):
            b = real[real["identity"].between(lo, hi)]
            r[f"{lo}-{hi}"] = ap.auc(b["truth"], b[col].fillna(-1))
        rows.append(r)
    return rows


def band_classifier(t, seed):
    sub = band_rows(t)
    y = sub["truth"].to_numpy()
    g = sub["meta_sample"].to_numpy()
    base = [c for c in ap.BASE if c in sub]
    print(f"  band {ap.BINS[0]}-{ap.BINS[-1]}, >= 3 fragments: {len(sub)} rows "
          f"({sub['class'].value_counts().to_dict()}); boosted classifier, samples grouped (5 folds):")
    variants = [("base", base), ("base + v17's three ancestry", base + ANC3), ("base + all five ancestry", base + ANC5)]
    if all(c in sub for c in FOREIGN):
        variants += [("base + foreign (full-reference scan)", base + FOREIGN),
                     ("base + five ancestry + foreign", base + ANC5 + FOREIGN)]
    for label, cols in variants:
        cols = [c for c in cols if c in sub]
        X = sub[cols].to_numpy(dtype=float)
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
    cols = [c for c in cols if c in sub]
    print("  medians of the errors against the right calls, rows with >= 3 fragments:")
    print(sub.groupby("error")[cols].median().T.round(4).to_string())
    e = sub[sub["error"].isin(["FP", "FN"])]
    rows = [{"feature": c, "AUC FN vs FP (FN high)": ap.auc((e["error"] == "FN").astype(int), e[c].fillna(-1))}
            for c in cols]
    print(pd.DataFrame(rows).round(3).to_string(index=False))


def scenario_errors(t):
    t = t.copy()
    y = t["truth"] == 1
    t["scenario"] = t["meta_scenario"].fillna("").astype(str).replace("", "design")
    t["FP"], t["FN"], t["TP"] = t["call"] & ~y, ~t["call"] & y, t["call"] & y
    g = t.groupby(["set", "scenario"])[["TP", "FP", "FN"]].sum()
    g["F1"] = 2 * g["TP"] / (2 * g["TP"] + g["FP"] + g["FN"]).clip(lower=1)
    print("  calls at the knob by collection and scenario (training rows: species held out; test: the final model):")
    print(g.round(4).to_string())
    fn = t[t["FN"]]
    print("  the FN by class and scenario (test and training):")
    print(pd.crosstab(fn["scenario"], fn["class"]).to_string())
    if "db_congeners_01" in t:
        print("  share of FN with a database congener within 0.01 (db_congeners_01 > 0), by scenario: " +
              ", ".join(f"{s} {v:.3f}" for s, v in fn.groupby("scenario")["db_congeners_01"].apply(lambda x: (x > 0).mean()).items()))
        tp = t[t["TP"]]
        print("  ... of the TP: " + ", ".join(f"{s} {v:.3f}" for s, v in
                                             tp.groupby("scenario")["db_congeners_01"].apply(lambda x: (x > 0).mean()).items()))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--build", required=True)
    p.add_argument("--v17")
    p.add_argument("--types", default="pe,se,pb,ont")
    p.add_argument("--seed", type=int, default=1)
    opts = p.parse_args(argv)
    pd.set_option("display.width", 250)
    for rt in opts.types.split(","):
        t = load(opts.build, rt)
        old = load(opts.v17, rt) if opts.v17 else None
        print(f"\n# {rt}: {len(t)} rows, {t['meta_sample'].nunique()} samples, knob {t['knob'].iloc[0]}\n", flush=True)
        print("## The error budget at the knob (v18)\n")
        ap.budget(t)
        print("\n## The ancestry features by class\n")
        by_class(t, "v18", ANC5 + FOREIGN + ["identity"])
        if old is not None:
            by_class(old, "v17", ANC3 + ["identity"])
        print("\n## AUC inside the identity band 0.95-0.985, present (high) against absent novel congener\n")
        rows = band_auc(t, "v18", ANC5 + FOREIGN + GAPS + OTHER)
        if old is not None:
            rows += band_auc(old, "v17", ANC3)
        print(pd.DataFrame(rows).round(3).to_string(index=False))
        print()
        band_classifier(t, opts.seed)
        print("\n## The errors on the features (v18)\n")
        errors(t, ANC5 + FOREIGN + ["identity", "lu_per_kb", "excess_scaled_median", "db_nearest_congener"])
        print("\n## The errors by scenario (v18)\n")
        scenario_errors(t)
    return 0


if __name__ == "__main__":
    sys.exit(main())
