#!/usr/bin/env python3
"""The soil false positives the models score highest, against the true calls they score as high (all features).

Rows: the soil and shallow-soil samples, hold-in scored with species held out and hold-out by the final model
(soil_errors.load). For each read type:
  1. the top false positives (score, fragments, identity, divergence, ...), and how often a taxon recurs as a false
     positive across the samples, and whether it is present (TP or FN) in other samples;
  2. false positives scored >= HIGH against true positives scored >= HIGH, matched on fragments (true positives
     drawn per fragment bin in the false positives' proportions): each feature's AUC (FP vs TP; 0.5 = no separation),
     sorted by how far from 0.5;
  3. a gradient-boosted classifier asked to separate those false positives from those true positives with every
     feature of the table (5 folds grouped by taxon): its AUC says whether the features hold the difference at all.

    python3 fp_top.py --build local/v13 --read-types pe,pb
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
from model_features import NON_FEATURE_COLUMNS  # noqa: E402

READ_TYPES = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
HIGH = 0.9
FRAG_BINS = [0, 2, 3, 5, 10, 30, 100, 1e12]
DIST_BINS = [-1, 0.0025, 0.005, 0.01, 0.015, 0.02, 0.03, 1]


def strain_kinds(rows):
    rep = pd.to_numeric(rows.meta_rep_genome, errors="coerce")
    ins = pd.to_numeric(rows.meta_insilico_strain, errors="coerce").fillna(0)
    kind = np.where(rep.eq(1), "representative", np.where(ins.eq(1), "in-silico strain", "real strain"))
    return pd.Series(kind).value_counts()


def load(build, rt):
    suffix = READ_TYPES[rt]
    name = "training_data.tsv" if rt == "pe" else f"training_data_{rt}.tsv"
    train = pd.read_csv(os.path.join(build, "training", name), sep="\t", low_memory=False,
                        dtype={"meta_scenario": str})
    pred = pd.read_csv(os.path.join(build, "model_logs", f"trained_model{suffix}.predictions.tsv.gz"), sep="\t",
                       usecols=["taxon", "p_species"])
    assert (pred.taxon.to_numpy() == train.taxon.to_numpy()).all()
    train["p"] = pred.p_species.to_numpy()
    train["set"] = "hold-in"
    test = pd.read_csv(os.path.join(build, "test", name), sep="\t", low_memory=False, dtype={"meta_scenario": str})
    test = test[test.meta_scenario.fillna("") != ""].copy()
    sp = pd.read_csv(os.path.join(build, "model_logs", f"trained_model{suffix}.scenario_predictions.tsv.gz"), sep="\t",
                     usecols=["taxon", "p"])
    assert (sp.taxon.to_numpy() == test.taxon.to_numpy()).all()
    test["p"] = sp.p.to_numpy()
    test["set"] = "hold-out"
    df = pd.concat([train, test], ignore_index=True)
    df = df[df.meta_scenario.isin(["soil", "soil_shallow"])].copy()
    with open(os.path.join(build, "model_logs", f"trained_model{suffix}.metrics.json")) as f:
        import json
        k = (json.load(f).get("depth_knobs") or {}).get("global_knob")
    df["knob"] = float(k) if k is not None else 0.5
    df["call"] = df.p >= df.knob
    df["outcome"] = np.select([df.truth.eq(1) & df.call, df.truth.eq(0) & df.call, df.truth.eq(1) & ~df.call],
                              ["TP", "FP", "FN"], "TN")
    return df


def feature_list(df):
    return [c for c in df.columns if c not in NON_FEATURE_COLUMNS and not c.startswith("meta_")
            and c not in ("p", "set", "knob", "call", "outcome", "domain") and pd.api.types.is_numeric_dtype(df[c])]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", required=True)
    ap.add_argument("--read-types", default="pe,pb")
    ap.add_argument("--top", type=int, default=25)
    opts = ap.parse_args()
    pd.set_option("display.width", 260)
    pd.set_option("display.max_columns", 40)
    pd.set_option("display.max_colwidth", 40)
    rng = np.random.default_rng(1)
    for rt in opts.read_types.split(","):
        df = load(opts.build, rt)
        feats = feature_list(df)
        print(f"\n######## {rt}: {len(df)} soil rows, knob {df.knob.iloc[0]}, {len(feats)} features")
        for sc in ("soil", "soil_shallow"):
            g = df[df.meta_scenario == sc]
            if not len(g):
                continue
            fp = g[g.outcome == "FP"]
            # recurrence: per taxon, samples where it is FP, and where it is present
            per = g.groupby("taxon").agg(fp_samples=("outcome", lambda s: (s == "FP").sum()),
                                         present_samples=("truth", "sum"), rows=("truth", "size"))
            fpt = per[per.fp_samples > 0]
            print(f"\n=== {rt} {sc}: {len(fp)} FP of {fpt.shape[0]} taxa; taxa FP in >= 2 samples: "
                  f"{(fpt.fp_samples >= 2).sum()}, >= 3: {(fpt.fp_samples >= 3).sum()}; FP taxa present in another "
                  f"sample: {(fpt.present_samples > 0).mean():.1%}; score quantiles of FP 50/75/90%: "
                  + "/".join(f"{q:.3f}" for q in fp.p.quantile([0.5, 0.75, 0.9])))
            cols = ["meta_sample", "taxon_name", "p", "fragments", "identity", "top_identity", "excess_scaled_median",
                    "hit_gene_fraction", "low_mapq_share", "congener_fit_share", "genus_share", "em_own_share",
                    "failed_candidate_rate", "mate_lost_share", "lu_per_kb", "adjacent_support", "cluster_genomes_log10"]
            cols = [c for c in cols if c in g]
            top = fp.sort_values("p", ascending=False).head(opts.top)[cols].copy()
            top["fp_samples"] = top.index.map(lambda i: per.fp_samples[g.taxon[i]])
            top["present_samples"] = top.index.map(lambda i: per.present_samples[g.taxon[i]])
            top["meta_sample"] = top.meta_sample.str.replace(r"^sc_", "", regex=True)
            print(top.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

            hi_fp = g[(g.outcome == "FP") & (g.p >= HIGH)]
            hi_tp = g[(g.outcome == "TP") & (g.p >= HIGH)]
            if len(hi_fp) < 20:
                continue
            fb = pd.cut(g.fragments, FRAG_BINS, right=False).astype(str)
            # matched on fragments only, and on fragments and the distance from the reference
            eb = pd.cut(g.excess_scaled_median.fillna(0), DIST_BINS).astype(str)
            for label, key in (("fragments", fb), ("fragments and distance", fb + " " + eb)):
                want = key[hi_fp.index].value_counts()
                picks, kept = [], []
                for b, n in want.items():
                    pool = hi_tp.index[key[hi_tp.index] == b]
                    if len(pool):
                        picks += list(rng.choice(pool, min(len(pool), 5 * n), replace=False))
                        kept += list(hi_fp.index[key[hi_fp.index] == b])
                tp_m, fp_m = g.loc[picks], g.loc[kept]
                both = pd.concat([fp_m, tp_m])
                yy = (both.outcome == "FP").to_numpy().astype(int)
                aucs = []
                for c in feats:
                    x = both[c].to_numpy(dtype=float)
                    if np.nanstd(x) == 0:
                        continue
                    a = roc_auc_score(yy, np.nan_to_num(x, nan=np.nanmedian(x)))
                    aucs.append((c, a, np.nanmedian(fp_m[c]), np.nanmedian(tp_m[c])))
                aucs.sort(key=lambda t: -abs(t[1] - 0.5))
                print(f"\nFP scored >= {HIGH} ({len(fp_m)}) vs TP scored >= {HIGH} matched on {label} ({len(tp_m)}; "
                      f"TP strain kinds: " + ", ".join(f"{k} {v}" for k, v in strain_kinds(tp_m).items()) +
                      "): the features that differ most (AUC: P(FP's value > TP's))")
                print(pd.DataFrame(aucs[:15], columns=["feature", "AUC", "FP median", "TP median"])
                      .to_string(index=False, float_format=lambda v: f"{v:.4f}"))
                X = both[feats].to_numpy(dtype=float)
                p = np.full(len(yy), np.nan)
                for tr, te in GroupKFold(5).split(X, yy, both.taxon.astype(str)):
                    m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                                       class_weight="balanced", random_state=1).fit(X[tr], yy[tr])
                    p[te] = m.predict_proba(X[te])[:, 1]
                print(f"a classifier of every feature separating them (5 folds by taxon): AUC {roc_auc_score(yy, p):.3f}")
            # among all called taxa: does a second model of every feature rank them better than the model's score?
            called = g[g.call]
            yy = (called.outcome == "FP").to_numpy().astype(int)
            X = called[feats].to_numpy(dtype=float)
            p2 = np.full(len(yy), np.nan)
            for tr, te in GroupKFold(5).split(X, yy, called.taxon.astype(str)):
                m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                                   class_weight="balanced", random_state=1).fit(X[tr], yy[tr])
                p2[te] = m.predict_proba(X[te])[:, 1]
            print(f"\namong the {len(called)} called taxa ({yy.sum()} FP): AUC for FP of the model's score "
                  f"{roc_auc_score(yy, -called.p):.3f}, of a second model of every feature fitted on soil's called "
                  f"taxa (5 folds by taxon) {roc_auc_score(yy, p2):.3f}")


if __name__ == "__main__":
    main()
