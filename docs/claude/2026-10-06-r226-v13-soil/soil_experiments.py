#!/usr/bin/env python3
"""Can another model, more features or the sample's context separate the soil scenarios' errors (r226 v13 tables)?

Each variant is fitted on a build's training table and scored two ways: with species held out (5 folds grouped by
taxon, as the trainer's p_species) on the training rows, and by the forest fitted on all of them on the test table's
rows (the design's test set and the scenarios' hold-out samples). Calls at a knob chosen as the trainer chooses its
global one (the threshold of the highest F1 with species held out over the rows the variant trains on, weighted as in
the fit, kept if it gains 0.002 over 0.5); the best F1 of each set at its own threshold beside it.

Variants (--variants):
  v13            the v13 trainer's forest: DEFAULT_FEATURE_SET, 64 trees, 512 leaves (pe, se) or 128 (pb, ont),
                 balanced classes, max_features sqrt, the scenarios' rows at weight 0.25
  weight 1       the scenarios' rows at weight 1
  soil only      fitted on the soil and soil_shallow hold-in rows alone (scored on soil only)
  +priors        DEFAULT_FEATURE_SET and the GTDB priors
  all            every feature column of the table
  big forest     256 trees, 4096 leaves
  gbm            histogram gradient boosting (500 rounds, rate 0.05, 63 leaves, balanced classes), same features
  +context       DEFAULT_FEATURE_SET and features of the taxon's sample and genus that protal could compute from the
                 profile (CONTEXT below)
  +context gbm   both
  forest 512     64 trees of 512 leaves for every read type (pb and ont: 4x their leaves)
  gbm large      gradient boosting with 127 leaves and 800 rounds
  +ref           DEFAULT_FEATURE_SET and the reference's k-mer uniqueness in the index (su_rate_ref, lu_rate_ref,
                 lsu_rate_ref: shares of the species' reference k-mers unique in the database; computed at build time)
  gbm depth rounded   gbm with sample_log_fragments rounded to 0.25 (a scenario sample's exact depth identifies it)
  gbm +ref depth rounded   both

CONTEXT (from the sample's rows; genus = the first word of the species name):
  sample_log_taxa         log10 of the taxa with reads in the sample
  sample_low_identity     the sample's fragments' share on taxa's low-identity bases (fragment-weighted
                          low_identity_share): how much of the sample is relatives the database lacks
  sample_identity         the fragment-weighted median identity of the sample's taxa with >= 10 fragments (read quality)
  identity_vs_sample      the taxon's identity minus sample_identity
  genus_log_ratio         log10((fragments + 1) / (genus_top_fragments + 1)): against its most abundant congener
  genus_log_taxa          log10 of the genus's taxa with reads in the sample
  genus_identity_gap      the taxon's identity minus the highest identity among its congeners with >= 3 fragments

    python3 soil_experiments.py --training local/v13/training --test local/v13/test --out results -t 6
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.model_selection import GroupKFold

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import machine_learning_cmdline as trainer  # noqa: E402
from model_features import DEFAULT_FEATURE_SET, PRIORS_FEATURES, feature_columns  # noqa: E402

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
LEAVES = {"pe": 512, "se": 512, "pb": 128, "ont": 128}
CONTEXT = ["sample_log_taxa", "sample_low_identity", "sample_identity", "identity_vs_sample", "genus_log_ratio",
           "genus_log_taxa", "genus_identity_gap"]
VARIANTS = ["v13", "weight 1", "soil only", "+priors", "all", "big forest", "gbm", "+context", "+context gbm"]
SOIL = ("soil", "soil_shallow")


def add_context(df):
    """The CONTEXT columns of a table's rows, from the other rows of their sample."""
    genus = df["taxon_name"].astype(str).str.replace("^s__", "", regex=True).str.split(" ").str[0]
    s = df["meta_sample"]
    frag = df["fragments"].astype(float)
    df["sample_log_taxa"] = np.log10(s.map(s.value_counts()).astype(float))
    w = frag * df["low_identity_share"].astype(float)
    df["sample_low_identity"] = w.groupby(s).transform("sum") / frag.groupby(s).transform("sum").clip(lower=1)
    deep = frag >= 10

    def weighted_median(g):
        g = g.sort_values("identity")
        c = g["fragments"].cumsum()
        return g["identity"].to_numpy()[np.searchsorted(c.to_numpy(), c.iloc[-1] / 2)]

    med = df[deep].groupby("meta_sample")[["identity", "fragments"]].apply(weighted_median)
    df["sample_identity"] = s.map(med).fillna(df["identity"].median()).astype(float)
    df["identity_vs_sample"] = df["identity"].astype(float) - df["sample_identity"]
    df["genus_log_ratio"] = np.log10((frag + 1) / (df["genus_top_fragments"].astype(float) + 1))
    key = s.astype(str) + "\t" + genus
    df["genus_log_taxa"] = np.log10(key.map(key.value_counts()).astype(float))
    # the highest identity among the taxon's congeners with >= 3 fragments (itself excluded)
    ident = df["identity"].astype(float).where(frag >= 3, -1.0)
    tmp = pd.DataFrame({"key": key, "ident": ident})
    top2 = tmp.groupby("key")["ident"].apply(lambda v: tuple(np.sort(v.to_numpy())[::-1][:2]) + (-1.0, -1.0))
    first = key.map(top2.map(lambda t: t[0]))
    second = key.map(top2.map(lambda t: t[1]))
    other = np.where(ident.to_numpy() >= first.to_numpy(), second.to_numpy(), first.to_numpy())
    df["genus_identity_gap"] = np.where(other >= 0, df["identity"].astype(float) - other, 0.05)
    return df


def model(variant, rt, seed, threads):
    if "gbm" in variant:
        rounds, leaves = (800, 127) if variant == "gbm large" else (500, 63)
        return HistGradientBoostingClassifier(max_iter=rounds, learning_rate=0.05, max_leaf_nodes=leaves,
                                              min_samples_leaf=20, l2_regularization=1.0, class_weight="balanced",
                                              early_stopping=False, random_state=seed)
    trees, leaves = {"big forest": (256, 4096), "forest 512": (64, 512)}.get(variant, (64, LEAVES[rt]))
    return RandomForestClassifier(n_estimators=trees, max_leaf_nodes=leaves, min_samples_leaf=1, max_features="sqrt",
                                  class_weight="balanced", random_state=seed, n_jobs=threads)


def f1_at(y, p, t, w=None):
    w = np.ones(len(y)) if w is None else w
    call = p >= t
    tp, fp, fn = (w * (call & (y == 1))).sum(), (w * (call & (y == 0))).sum(), (w * (~call & (y == 1))).sum()
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0


def choose_knob(y, p, w):
    grid = np.round(np.arange(0.30, 0.851, 0.01), 2)
    f = np.array([f1_at(y, p, t, w) for t in grid])
    base = f1_at(y, p, 0.5, w)
    return float(grid[f.argmax()]) if f.max() >= base + 0.002 else 0.5


def best(y, p):
    grid = np.round(np.arange(0.05, 0.96, 0.01), 2)
    f = np.array([f1_at(y, p, t) for t in grid])
    return float(f.max()), float(grid[f.argmax()])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--training", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("-t", "--threads", type=int, default=6)
    opts = ap.parse_args()
    warnings.simplefilter("ignore", pd.errors.PerformanceWarning)
    os.makedirs(opts.out, exist_ok=True)
    rows_out = []
    for rt in opts.read_types.split(","):
        began = time.time()
        train = add_context(trainer.load_table(os.path.join(opts.training, TABLES[rt])))
        test = add_context(trainer.load_table(os.path.join(opts.test, TABLES[rt])))
        y, yt = train["truth"].to_numpy(), test["truth"].to_numpy()
        sc, sct = trainer.scenario_of(train), trainer.scenario_of(test)
        base = feature_columns(train.columns, DEFAULT_FEATURE_SET)
        everything = [c for c in feature_columns(train.columns, "all")
                      if c not in CONTEXT and pd.api.types.is_numeric_dtype(train[c]) and c != "domain"]
        splits = list(GroupKFold(opts.folds).split(train, y, train["taxon"].astype(str)))
        sets = [("design, species held out", "train", sc == ""), ("soil hold-in, species held out", "train", sc == "soil"),
                ("soil_shallow hold-in, species held out", "train", sc == "soil_shallow"),
                ("design test", "test", sct == ""), ("gut hold-out", "test", sct == "gut"),
                ("soil hold-out", "test", sct == "soil"), ("soil_shallow hold-out", "test", sct == "soil_shallow")]
        for variant in opts.variants.split(","):
            ref = ["su_rate_ref", "lu_rate_ref", "lsu_rate_ref"]
            cols = {"+priors": base + [c for c in PRIORS_FEATURES if c not in base], "all": everything,
                    "+context": base + CONTEXT, "+context gbm": base + CONTEXT, "+ref": base + ref,
                    "gbm +ref depth rounded": base + ref}.get(variant, base)
            X, Xt = train[cols].to_numpy(dtype=np.float64), test[cols].to_numpy(dtype=np.float64)
            if "depth rounded" in variant and "sample_log_fragments" in cols:
                j = cols.index("sample_log_fragments")
                X[:, j], Xt[:, j] = np.round(X[:, j] * 4) / 4, np.round(Xt[:, j] * 4) / 4
            fit_rows = np.isin(sc, SOIL) if variant == "soil only" else np.ones(len(y), dtype=bool)
            w_all = np.where(sc == "", 1.0, 1.0 if variant in ("weight 1", "soil only") else 0.25)
            p = np.full(len(y), np.nan)
            for tr, te in splits:
                tr = tr[fit_rows[tr]]
                te = te[fit_rows[te]]
                m = model(variant, rt, opts.seed, opts.threads).fit(X[tr], y[tr], sample_weight=w_all[tr])
                p[te] = m.predict_proba(X[te])[:, 1]
            m = model(variant, rt, opts.seed, opts.threads).fit(X[fit_rows], y[fit_rows], sample_weight=w_all[fit_rows])
            pt = m.predict_proba(Xt)[:, 1]
            knob = choose_knob(y[fit_rows], p[fit_rows], w_all[fit_rows])
            for label, where, mask in sets:
                yy, pp = (y[mask], p[mask]) if where == "train" else (yt[mask], pt[mask])
                ok = ~np.isnan(pp)
                if not ok.any():
                    continue
                yy, pp = yy[ok], pp[ok]
                call = pp >= knob
                bf, bt = best(yy, pp)
                rows_out.append({"read type": rt, "variant": variant, "set": label, "knob": knob,
                                 "F1": f1_at(yy, pp, knob), "FP": int((call & (yy == 0)).sum()),
                                 "FN": int((~call & (yy == 1)).sum()), "F1 at 0.5": f1_at(yy, pp, 0.5),
                                 "best F1": bf, "at": bt})
            print(f"{rt} {variant}: knob {knob}, {time.time() - began:.0f} s in; " +
                  ", ".join(f"{r['set'].split(',')[0]} {r['F1']:.4f}" for r in rows_out if r["read type"] == rt and r["variant"] == variant),
                  flush=True)
            pd.DataFrame(rows_out).to_csv(os.path.join(opts.out, "runs.tsv"), sep="\t", index=False, float_format="%.5g")
    out = pd.DataFrame(rows_out)
    wide = out.pivot_table(index=["read type", "set"], columns="variant", values="F1", sort=False)
    wide = wide[[v for v in opts.variants.split(",") if v in wide.columns]]
    pd.set_option("display.width", 300)
    print(wide.to_string(float_format=lambda v: f"{v:.4f}"))
    wide.to_csv(os.path.join(opts.out, "summary.tsv"), sep="\t", float_format="%.5g")


if __name__ == "__main__":
    main()
