"""Calls of species the database lacks reported as "an unknown species of genus G" instead of a congener.

On the taxa a model calls at its knob curve (as protal does), a second forest scores whether the taxon's reads are a
congener's the database lacks: target 1 for an absent taxon of a genus with a species in the sample that the database lacks
(meta_novel_congener), 0 otherwise. Called taxa scoring at or above a threshold are converted: no species call, one
"unknown species of genus G" call per sample and genus instead. The threshold is chosen on the training rows by
out-of-fold scores (5 folds grouped by species, so a fold never saw the taxon), the forest then fitted on all training rows
and applied once to the independent test set.

Species level: a converted present taxon becomes a false negative, a converted absent one is no longer a false positive.
Genus level: an "unknown species of G" call is right if the sample holds a species of G the database lacks.

Usage: genus_fallback.py RUN_DIR [SUFFIX ...], e.g. local/v5 "" _se _pb _ont (RUN_DIR holds training/, test/, the
model's metrics.json, predictions and test predictions)."""
import json
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "scripts"))
from model_features import ADJACENCY_FEATURES, NORMALIZED_FEATURES, RELATIVE_FEATURES  # noqa: E402

FEATURES = NORMALIZED_FEATURES + ADJACENCY_FEATURES + RELATIVE_FEATURES + ["p", "sample_log_fragments"]
THRESHOLDS = np.round(np.arange(0.05, 1.0, 0.05), 2)


def f1(tp, fp, fn):
    return 2 * tp / max(1, 2 * tp + fp + fn)


def load(run, part, suffix, metrics):
    """The rows of training/ or test/ with the model's score p (species held out for training) and the knob curve's call."""
    name = "predictions" if part == "training" else "test_predictions"
    pred = pd.read_csv(os.path.join(run, f"trained_model{suffix}.{name}.tsv.gz"), sep="\t", low_memory=False)
    if part == "training":
        pred = pred.rename(columns={"p_species": "p"})
    pred = pred[["meta_sample", "taxon", "p"]]
    table = pd.read_csv(os.path.join(run, part, f"training_data{suffix}.tsv"), sep="\t", low_memory=False)
    d = table.merge(pred, on=["meta_sample", "taxon"], how="inner", validate="one_to_one")
    assert len(d) == len(table), f"{part}: {len(table) - len(d)} rows without a score"
    lineage = pd.read_csv(os.path.join(run, f"trained_model{suffix}.{name}.tsv.gz"), sep="\t", low_memory=False,
                          usecols=["meta_sample", "taxon", "meta_lineage_genus"])
    d = d.merge(lineage, on=["meta_sample", "taxon"], how="left").copy()
    sample_fragments = d.groupby("meta_sample")["fragments"].transform("sum")
    d["sample_log_fragments"] = np.log10(np.maximum(sample_fragments, 1))
    curve = np.array(metrics["depth_knobs"]["curve"], dtype=float)
    d["called"] = d["p"] >= np.interp(d["sample_log_fragments"], curve[:, 0], curve[:, 1])
    d["present"] = d["truth"] == 1
    d["novel_genus"] = d["meta_novel_congener"] == 1
    d["target"] = ~d["present"] & d["novel_genus"]
    return d


def outcome(d, convert):
    """Species-level counts with `convert` (rows) turned into genus calls, and the genus calls' correctness."""
    call = d["called"] & ~convert
    tp, fp = int((call & d.present).sum()), int((call & ~d.present).sum())
    fn = int((~call & d.present).sum())
    conv = d[convert]
    genus_calls = conv.groupby(["meta_sample", "meta_lineage_genus"])["novel_genus"].first()
    return {"F1": f1(tp, fp, fn), "TP": tp, "FP": fp, "FN": fn, "converted": len(conv),
            "converted absent, genus right": int((~conv.present & conv.novel_genus).sum()),
            "converted absent, genus wrong": int((~conv.present & ~conv.novel_genus).sum()),
            "converted present (species lost)": int(conv.present.sum()),
            "  of them strains": int((conv.present & (conv.meta_rep_genome == 0)).sum()),
            "  of them with a missing congener too": int((conv.present & conv.novel_genus).sum()),
            "genus calls": len(genus_calls), "genus calls right": int(genus_calls.sum())}


def print_separation(tr, oof, suffix):
    """How well the out-of-fold scores separate the target among the called training rows, and the groups' medians."""
    t = suffix.strip("_") or "pe"
    print(f"## {t}: out of fold on the called training rows: AUC {roc_auc_score(tr.target, oof):.3f}, "
          f"AP {average_precision_score(tr.target, oof):.3f} (base rate {tr.target.mean():.3f})")
    groups = {"absent, genus has a missing species (target)": tr.target,
              "absent, other": ~tr.present & ~tr.novel_genus,
              "present, representative": tr.present & (tr.meta_rep_genome == 1),
              "present, strain": tr.present & (tr.meta_rep_genome == 0),
              "present, strain, 1-10 fragments": tr.present & (tr.meta_rep_genome == 0) & (tr.fragments <= 10)}
    cols = ["fragments", "identity", "top_identity", "excess_median", "genus_spill", "em_own_share", "p"]
    medians = pd.DataFrame({k: tr.loc[m, cols].median() for k, m in groups.items()}).T
    medians["rows"] = [int(m.sum()) for m in groups.values()]
    print(medians.to_string(float_format=lambda v: f"{v:.3f}"))


def forest(seed=1):
    return RandomForestClassifier(n_estimators=300, min_samples_leaf=3, max_features="sqrt",
                                  class_weight="balanced_subsample", n_jobs=-1, random_state=seed)


def run_type(run, suffix):
    metrics = json.load(open(os.path.join(run, f"trained_model{suffix}.metrics.json")))
    train, test = load(run, "training", suffix, metrics), load(run, "test", suffix, metrics)
    features = [c for c in FEATURES if c in train.columns]
    tr = train[train.called].reset_index(drop=True)
    # Out-of-fold scores on the called training rows, folds grouped by species.
    oof = np.zeros(len(tr))
    for fit_idx, score_idx in GroupKFold(n_splits=5).split(tr, tr.target, groups=tr.taxon):
        m = forest().fit(tr.loc[fit_idx, features], tr.loc[fit_idx, "target"])
        oof[score_idx] = m.predict_proba(tr.loc[score_idx, features])[:, 1]
    scores_train = pd.Series(0.0, index=train.index)
    scores_train[train.called] = oof
    print_separation(tr, oof, suffix)
    best = max(THRESHOLDS, key=lambda t: outcome(train, train.called & (scores_train >= t))["F1"])
    m = forest().fit(tr[features], tr.target)
    scores_test = pd.Series(0.0, index=test.index)
    scores_test[test.called] = m.predict_proba(test.loc[test.called, features])[:, 1]

    t = suffix.strip("_") or "pe"
    print(f"## {t}: {run}, {len(tr)} called training rows ({int(tr.target.sum())} absent of a genus with a missing species), "
          f"{int(test.called.sum())} called test rows ({int((test.called & test.target).sum())}), threshold {best} "
          f"(chosen on the training rows, out of fold)")
    keys = ["F1", "TP", "FP", "FN", "converted", "converted absent, genus right", "converted absent, genus wrong",
            "converted present (species lost)", "  of them strains", "  of them with a missing congener too",
            "genus calls", "genus calls right"]
    rows = {"training, no conversion": outcome(train, train.called & False),
            f"training, out of fold at {best}": outcome(train, train.called & (scores_train >= best)),
            "test, no conversion": outcome(test, test.called & False),
            f"test at {best}": outcome(test, test.called & (scores_test >= best))}
    print(pd.DataFrame(rows).loc[keys].to_string(float_format=lambda v: f"{v:.4f}"))
    n_samples = test.meta_sample.nunique()
    conv = outcome(test, test.called & (scores_test >= best))
    print(f"test: {conv['converted'] / n_samples:.2f} species converted per sample ({n_samples} samples), "
          f"{conv['genus calls'] / n_samples:.2f} genus calls; species F1 {rows['test, no conversion']['F1']:.4f} -> "
          f"{conv['F1']:.4f}")
    print("test, by threshold (the chosen one is marked; the others are for the trade-off, not choices):")
    sweep = {}
    for th in THRESHOLDS:
        o = outcome(test, test.called & (scores_test >= th))
        sweep[f"{th}{' *' if th == best else ''}"] = {k: o[k] for k in ("F1", "converted", "converted absent, genus right",
                                                                        "converted absent, genus wrong",
                                                                        "converted present (species lost)",
                                                                        "genus calls right", "genus calls")}
    print(pd.DataFrame(sweep).T.to_string(float_format=lambda v: f"{v:.4f}"))
    by_depth = test[test.called].assign(conv=(scores_test >= best)[test.called].to_numpy())
    print("test at the chosen threshold, by depth: converted (absent, genus right / present lost)")
    for depth, g in by_depth.groupby("meta_read_pairs"):
        c = g[g.conv]
        print(f"  {depth:>10}: {len(c):3} ({int((~c.present & c.novel_genus).sum())} / {int(c.present.sum())}) of {len(g)} "
              f"called, {int((~g.present).sum())} of them false")
    imp = pd.Series(m.feature_importances_, index=features).sort_values(ascending=False).head(8)
    print("importance: " + ", ".join(f"{k} {v:.3f}" for k, v in imp.items()))
    print()


if __name__ == "__main__":
    for s in (sys.argv[2:] or [""]):
        run_type(sys.argv[1], s)
