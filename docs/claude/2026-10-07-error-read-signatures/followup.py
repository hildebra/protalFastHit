#!/usr/bin/env python3
"""Follow-ups to signatures.py: the read-level model with at most 20 fragments per taxon and without the reads'
sequence composition (GC and the like name the species, not the error), with permutation importances; and whether
a taxon's per-gene divergence fits the within-species or the between-species pattern of gene_congeners.tsv.

    python3 followup.py --records ~/v15/err_analysis --build local/v15 > followup.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold

import signatures

COMPOSITION = ["gc", "top_trinucleotide", "homopolymer"]
CAP = 20


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--types", default="pe,se,pb,ont")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


def capped(f, seed):
    """At most CAP fragments per sample, taxon and group."""
    return (f.sample(frac=1.0, random_state=seed).groupby(["sample", "taxon", "group"], sort=False).head(CAP))


def read_model(f, seed):
    c = capped(f[f["counted"] & f["group"].isin(["wrong (FP)", "own (FN)", "own (TP)"])], seed)
    y = (c["group"] == "wrong (FP)").astype(int).to_numpy()
    cols = [x for x in signatures.FEATURES if x in c and x not in COMPOSITION and not c[x].isna().all()
            and c[x].astype(float).nunique() > 1]
    X = c[cols].astype(float)
    groups = c["sample"].to_numpy()
    s = np.zeros(len(c))
    for tr, te in GroupKFold(n_splits=5).split(X, y, groups):
        m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, random_state=seed)
        m.fit(X.iloc[tr], y[tr])
        s[te] = m.predict_proba(X.iloc[te])[:, 1]
    print(f"capped at {CAP} fragments per taxon: {len(c)} fragments ({y.sum()} FP, "
          f"{(c['group'] == 'own (FN)').sum()} own FN, {(c['group'] == 'own (TP)').sum()} own TP)")
    print(f"boosted, no composition, sample-grouped 5-fold: AUC {roc_auc_score(y, s):.3f}; "
          f"FP vs own FN {roc_auc_score(y[c['group'] != 'own (TP)'], s[c['group'] != 'own (TP)']):.3f}; "
          f"FP vs own TP {roc_auc_score(y[c['group'] != 'own (FN)'], s[c['group'] != 'own (FN)']):.3f}")
    rows = []
    for col in cols:
        rows.append({"feature": col, "AUC": signatures.auc(y, X[col]), "FP": X.loc[y == 1, col].mean(),
                     "own FN": X.loc[c["group"].to_numpy() == "own (FN)", col].mean(),
                     "own TP": X.loc[c["group"].to_numpy() == "own (TP)", col].mean()})
    m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, random_state=seed).fit(X, y)
    sub = np.random.RandomState(seed).choice(len(X), min(len(X), 20000), replace=False)
    imp = permutation_importance(m, X.iloc[sub], y[sub], scoring="roc_auc", n_repeats=3, random_state=seed)
    out = pd.DataFrame(rows)
    out["importance"] = imp.importances_mean
    out = out.sort_values("importance", ascending=False)
    signatures.show(out.head(14).round(3), index=False)
    # The composition's AUC, for the record
    for col in COMPOSITION:
        if col in c:
            print(f"  ({col}: AUC {signatures.auc(y, c[col]):.3f}, left out)")


def gene_pattern(f, genes):
    """Per taxon with >= 8 counted fragments on >= 4 genes: the slope of fragment divergence on the genes' within- and
    between-species factors, and which of the two explains it better (R^2)."""
    c = f[f["counted"] & f["role"].isin(["FP", "FN"]) & f["gene"].notna()].copy()
    c["div"] = 1 - c["identity"]
    w = genes["within_factor"].to_dict()
    b = genes["between_factor"].to_dict()
    c["fw"], c["fb"] = c["gene"].map(w), c["gene"].map(b)
    print(f"across the 168 genes: Spearman(within_factor, between_factor) = "
          f"{genes['within_factor'].corr(genes['between_factor'], method='spearman'):.2f}")
    rows = []
    for (sample, taxon, role), g in c.groupby(["sample", "taxon", "role"]):
        if len(g) < 8 or g["gene"].nunique() < 4:
            continue
        gm = g.groupby("gene").agg(div=("div", "mean"), fw=("fw", "first"), fb=("fb", "first"), n=("div", "size"))
        r2 = {}
        for col in ("fw", "fb"):
            x, y, wt = gm[col].to_numpy(), gm["div"].to_numpy(), gm["n"].to_numpy()
            slope = (wt * x * y).sum() / max(1e-12, (wt * x * x).sum())  # through the origin: div = slope * factor
            resid = y - slope * x
            r2[col] = 1 - (wt * resid ** 2).sum() / max(1e-12, (wt * (y - np.average(y, weights=wt)) ** 2).sum())
        rows.append({"role": role, "fragments": len(g), "genes": len(gm), "div": g["div"].mean(),
                     "r2_within": r2["fw"], "r2_between": r2["fb"], "between_better": r2["fb"] - r2["fw"],
                     "cv": gm["div"].std() / max(1e-9, gm["div"].mean())})
    d = pd.DataFrame(rows)
    if d.empty:
        return
    d["truth"] = (d["role"] == "FN").astype(int)
    print(f"{len(d)} taxa ({(d['role'] == 'FP').sum()} FP, {(d['role'] == 'FN').sum()} FN)")
    signatures.show(d.groupby("role")[["fragments", "genes", "div", "r2_within", "r2_between", "between_better", "cv"]]
                    .median().round(4))
    for col in ("between_better", "cv", "div"):
        print(f"  AUC of {col} for FN (present) against FP: {signatures.auc(d['truth'], d[col]):.3f}")


def fp_other(f):
    """The FP fragments from beyond the genus: identity and source."""
    w = f[f["counted"] & (f["group"] == "wrong (FP)")]
    for rel in ("genus", "family", "order", "other"):
        g = w[w["relation"] == rel]
        if len(g):
            print(f"  {rel:7s} {len(g):7d} fragments, {g.groupby(['sample', 'taxon']).ngroups:5d} taxa: identity median "
                  f"{g['identity'].median():.3f}, >= 0.99: {(g['identity'] >= 0.99).mean():.3f}; source held out "
                  f"{g['source_held_out'].mean():.3f}; MAPQ median {g['mapq'].median():.0f}")


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    genes = pd.read_csv(os.path.join(opts.build, "model_logs", "gene_congeners.tsv"), sep="\t").set_index("geneid")
    genes.index = genes.index.astype(float)
    for rt in opts.types.split(","):
        f = pd.read_pickle(os.path.join(opts.records, f"{rt}.fragments.pkl.gz"))
        f["group"] = signatures.group_of(f)
        print(f"\n# {rt}\n\n## Read level, capped and without composition\n", flush=True)
        read_model(f, opts.seed)
        print("\n## FP fragments by their source's relation\n")
        fp_other(f)
        print("\n## Per-gene divergence: within- or between-species pattern\n")
        gene_pattern(f, genes)
    return 0


if __name__ == "__main__":
    sys.exit(main())
