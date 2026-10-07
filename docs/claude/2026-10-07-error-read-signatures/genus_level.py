#!/usr/bin/env python3
"""Genus-level calls for novel congeners, on the v15 build's own calls: what the false positives at the knob are, the
ceiling if every false positive that is a congener of a species the database lacks were reported as "an unknown species
of genus G" (a right genus-level call), and what a classifier on protal's features would convert, its threshold chosen
out of fold on the training rows (samples grouped) and applied once to the test rows.

    python3 genus_level.py --records ~/v15/err_analysis --build local/v15 > genus_level.txt
"""
import argparse
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

import calls

THRESHOLDS = np.round(np.arange(0.2, 0.96, 0.05), 2)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--types", default="pe,se,pb,ont")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


def budget(label, c):
    """The calls of one set: errors by class, species F1 and the ceiling of genus-level reporting."""
    called = c["call"]
    y = c["truth"] == 1
    tp, fp, fn = int((called & y).sum()), int((called & ~y).sum()), int((~called & y).sum())
    fp_by = c.loc[called & ~y, "class"].value_counts()
    novel = int(fp_by.get("absent novel congener", 0))
    fn_by = c.loc[~called & y, "class"].value_counts()
    print(f"{label}: {c['meta_sample'].nunique()} samples, {len(c)} rows, {int(y.sum())} present; called {int(called.sum())}: "
          f"TP {tp}, FP {fp}, FN {fn}, F1 {calls.f1(tp, fp, fn):.4f}")
    print(f"  FP by class: {fp_by.to_dict()}; FN by class: {fn_by.to_dict()}")
    print(f"  ceiling, the {novel} novel-congener FP reported at genus level: F1 {calls.f1(tp, fp - novel, fn):.4f} "
          f"(FP dropped), {calls.f1(tp + novel, fp - novel, fn):.4f} (counted as right genus calls)")
    # The genus of such a false positive is right by definition; how often a false positive's genus is right at all:
    # the FP of the other absent classes have their closest sample species in another genus or in the database.
    print(f"  genus right if every FP were converted: {novel / max(1, fp):.3f}")


def convert(c, scores, thr, fn_uncalled=0):
    """Species-level errors when called taxa scoring >= thr become genus calls: (converted present, converted absent
    novel congener, converted absent other, F1 with the converted absent dropped, F1 with them counted right).
    c: called rows; fn_uncalled: the set's present taxa not called, which stay false negatives."""
    called = c["call"].to_numpy()
    y = (c["truth"] == 1).to_numpy()
    cls = c["class"].to_numpy()
    conv = called & (scores >= thr)
    tp = int((called & ~conv & y).sum())
    fp = int((called & ~conv & ~y).sum())
    fn = int(((~called | conv) & y).sum()) + fn_uncalled
    right = int((conv & (cls == "absent novel congener")).sum())
    return (int((conv & y).sum()), right, int((conv & ~y).sum()) - right, calls.f1(tp, fp, fn),
            calls.f1(tp + right, fp, fn))


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    for rt in opts.types.split(","):
        t, features = calls.load(opts.build, rt)
        knob = calls.knob_of(opts.records, rt)
        t["call"] = t["p"] >= knob
        t["class"] = calls.classes(t)
        t["set"] = t["meta_sample"].str.split(":").str[0]
        print(f"\n# {rt}: knob {knob}, {len(features)} features\n", flush=True)
        for label, c in (("training rows, species held out", t[t["set"] == "training"]),
                         ("test rows (design test set + scenario hold-out samples where kept)", t[t["set"] == "test"])):
            budget(label, c)

        print("\n## A classifier for 'absent novel congener' among the called taxa\n", flush=True)
        c = t[t["call"]].reset_index(drop=True)
        uncalled = {s: int((~t["call"] & (t["truth"] == 1) & (t["set"] == s)).sum()) for s in ("training", "test")}
        y = (c["class"] == "absent novel congener").astype(int).to_numpy()
        X = c[features].to_numpy(dtype=float)
        X = np.where(np.isfinite(X), X, np.nan)
        train = (c["set"] == "training").to_numpy()
        groups = c.loc[train, "meta_sample"].to_numpy()
        oof = np.full(train.sum(), np.nan)
        for tr, te in GroupKFold(n_splits=5).split(X[train], y[train], groups):
            m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.1, max_leaf_nodes=31, random_state=opts.seed)
            m.fit(X[train][tr], y[train][tr])
            oof[te] = m.predict_proba(X[train][te])[:, 1]
        ytr = y[train]
        print(f"called training rows {train.sum()}, targets {int(ytr.sum())}: out-of-fold AUC "
              f"{roc_auc_score(ytr, oof):.3f}, AP {average_precision_score(ytr, oof):.3f}")
        rows = []
        ctr = c[train].reset_index(drop=True)
        base_drop, base_right = convert(ctr, oof, 2.0, uncalled["training"])[3:]
        for thr in THRESHOLDS:
            pres, right, wrong, f_drop, f_right = convert(ctr, oof, thr, uncalled["training"])
            rows.append({"threshold": thr, "converted present": pres, "absent, genus right": right,
                         "absent, genus wrong": wrong, "F1 (FP dropped)": round(f_drop, 4),
                         "F1 (genus right counted)": round(f_right, 4)})
        table = pd.DataFrame(rows)
        print(f"training rows out of fold (no conversion: F1 {base_drop:.4f}):")
        print(table.to_string(index=False))
        best = table.loc[table["F1 (FP dropped)"].idxmax(), "threshold"]
        best_r = table.loc[table["F1 (genus right counted)"].idxmax(), "threshold"]
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.1, max_leaf_nodes=31, random_state=opts.seed)
        m.fit(X[train], ytr)
        test = ~train
        if test.sum():
            s = m.predict_proba(X[test])[:, 1]
            cte = c[test].reset_index(drop=True)
            yte = y[test]
            print(f"\ntest rows: called {test.sum()}, targets {int(yte.sum())}, AUC {roc_auc_score(yte, s):.3f}, "
                  f"AP {average_precision_score(yte, s):.3f}; no conversion F1 "
                  f"{convert(cte, s, 2.0, uncalled['test'])[3]:.4f}")
            for label, thr in (("threshold chosen on training rows, FP dropped", best),
                               ("threshold chosen on training rows, genus right counted", best_r)):
                pres, right, wrong, f_drop, f_right = convert(cte, s, thr, uncalled["test"])
                print(f"  {label} ({thr}): converted present {pres}, absent genus right {right}, wrong {wrong}; "
                      f"F1 {f_drop:.4f} (FP dropped), {f_right:.4f} (genus right counted)")
            rows = [dict(threshold=thr, **dict(zip(["converted present", "absent, genus right", "absent, genus wrong",
                                                    "F1 (FP dropped)", "F1 (genus right counted)"],
                                                   convert(cte, s, thr, uncalled["test"])))) for thr in THRESHOLDS]
            print("  the trade-off on the test rows (not a choice):")
            print(pd.DataFrame(rows).round(4).to_string(index=False))
        imp = pd.Series(m.feature_importances_ if hasattr(m, "feature_importances_") else np.zeros(len(features)),
                        index=features)
        if imp.sum() == 0:
            from sklearn.inspection import permutation_importance
            sub = np.random.RandomState(opts.seed).choice(train.sum(), min(4000, train.sum()), replace=False)
            r = permutation_importance(m, X[train][sub], ytr[sub], scoring="roc_auc", n_repeats=3,
                                       random_state=opts.seed)
            imp = pd.Series(r.importances_mean, index=features)
        print("\n  top features (permutation importance on training rows):")
        print(imp.sort_values(ascending=False).head(12).round(4).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
