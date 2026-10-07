#!/usr/bin/env python3
"""What tells the reads of false positives and false negatives apart, from error_read_features.py's tables.

Reads OUT/<set>/<point>/<sample>.taxa.tsv.gz (one line per row: TP, FP, FN, TN) and .records.tsv.gz (the records of
the FP and FN fragments) and prints:
  1. the rows by class, set and scenario
  2. FP anatomy: where the fragments counted for FP taxa come from (own species, same genus, family, other; held out)
     against the TP taxa's
  3. FN anatomy: where the reads of the FN taxa's own species went (counted here, low MAPQ here, elsewhere, unaligned)
  4. each taxon feature's separation of truth where the model decides: AUC among called rows (TP vs FP), among uncalled
     rows (FN vs TN), and within bins of the counted fragments and of the model's score p
  5. whether the features add to the model's score: sample-grouped cross-validation of truth ~ p against
     truth ~ p + features (gradient boosting), log-loss, AUC, best F1 and its FP/FN, permutation importances
  6. read level: among the counted records of the FP and FN fragments, what tells an alignment to another species than
     the read's own (the FP taxa's reads) from an alignment to its own (the FN taxa's)

    python3 analyse_error_reads.py --dir local/v15/error_read_features/pe > analysis_pe.txt

Needs pandas, numpy and scikit-learn.
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss, roc_auc_score
from sklearn.model_selection import GroupKFold

# Taxon features from the reads (taxa table), and ratios made of its counts.
SHARES = ["identity_mean", "identity_min", "identity_sd", "mapq_mean", "mapq_max_share", "zu_mean", "zt_mean",
          "zu0_share", "zt0_share", "za_share", "za_tie_share", "za_close_share", "za_called_share", "zf_share",
          "zf_mean", "zf_called_share", "clip_share", "edge_share", "both_mates_share", "split_genes_share",
          "other_taxon_share", "proper_share", "mate_unmapped_share", "zr_share", "secondary_share", "aln_mean", "genes", "fragments_per_gene_max"]
RATIOS = {
    "lowmapq_share": ("frags_lowmapq", "frags_best"),          # best on the taxon but below MAPQ 4
    "any_not_best_share": ("any_not_best", "frags_any"),       # a record on the taxon, the best elsewhere
    "failed_per_counted": ("seeded_failed", "frags_counted1"),  # reads that seeded on it and failed, per counted one
    "own_elsewhere_share": ("own_elsewhere", "own_frags"),      # oracle: its own species' reads that went elsewhere
}
RECORD_FEATURES = ["identity", "mapq", "zu", "zt", "zu_per_100", "za_n", "za_min_filled", "za_ties", "zf_n",
                   "clipped", "at_edge", "proper", "abs_tlen", "frag_aligned_taxa", "aln_query", "gap_opens",
                   "mismatches", "za_called", "zf_called"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--dir", required=True, help="error_read_features.py's output folder (one read type)")
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


def load(folder, kind):
    paths = sorted(glob.glob(os.path.join(folder, "*", "*", f"*.{kind}.tsv.gz")))
    if not paths:
        sys.exit(f"no *.{kind}.tsv.gz under {folder}")
    return pd.concat((pd.read_csv(p, sep="\t", low_memory=False) for p in paths), ignore_index=True)


def auc(y, x):
    """AUC of x for y (0.5: none; below 0.5: lower in the positives), or NaN."""
    ok = ~pd.isna(x)
    y, x = np.asarray(y)[ok], np.asarray(x, dtype=float)[ok]
    if len(set(y)) < 2:
        return float("nan")
    return roc_auc_score(y, x)


def binned_auc(df, feature, by):
    """The AUC of a feature for truth within bins of `by`, weighted by the bins' smaller class."""
    total, weight = 0.0, 0
    for _, g in df.groupby(by, observed=True):
        a = auc(g["truth"], g[feature])
        if np.isnan(a):
            continue
        w = min(g["truth"].sum(), len(g) - g["truth"].sum())
        total, weight = total + a * w, weight + w
    return total / weight if weight else float("nan")


def section(title):
    print(f"\n## {title}\n", flush=True)


def taxa_frame(folder):
    t = load(folder, "taxa")
    t["truth"] = t["truth"].astype(str).str.lower().isin(["1", "true"]).astype(int)
    t["p"] = pd.to_numeric(t["p"], errors="coerce")
    t["any_not_best"] = (t["frags_any"] - t["frags_best"]).clip(lower=0)
    t["frags_counted1"] = t["frags_counted"].clip(lower=1)
    t["own_elsewhere"] = t["own_counted_elsewhere"] + t["own_lowmapq_elsewhere"]
    for name, (num, den) in RATIOS.items():
        t[name] = np.where(t[den] > 0, t[num] / t[den].where(t[den] > 0, 1), np.nan)
    t["frag_bin"] = pd.cut(t["frags_counted"], [-1, 0, 1, 2, 5, 20, 1e12], labels=["0", "1", "2", "3-5", "6-20", ">20"])
    t["p_bin"] = pd.cut(t["p"], [-0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 1.01])
    return t


def overview(t):
    section("1. Rows by class")
    print(pd.crosstab([t["set"], t["scenario"]], t["class"]).to_string())


def fp_anatomy(t):
    section("2. Where the counted fragments of FP and TP taxa come from (shares of all their counted fragments)")
    cols = ["src_own", "src_genus", "src_family", "src_other", "src_host", "src_unknown", "src_held_out",
            "src_strain", "src_insilico"]
    rows = []
    for (cls, scenario), g in t[t["class"].isin(["FP", "TP"])].groupby(["class", "scenario"]):
        n = g["frags_counted"].sum()
        rows.append(dict({"class": cls, "scenario": scenario, "taxa": len(g), "fragments": n,
                          "median fragments": g["frags_counted"].median()},
                         **{c[4:]: round(g[c].sum() / n, 3) if n else np.nan for c in cols}))
    print(pd.DataFrame(rows).to_string(index=False))
    fp = t[t["class"] == "FP"].copy()
    if len(fp):
        fp["main"] = fp[["src_own", "src_genus", "src_family", "src_other", "src_host", "src_unknown"]].idxmax(axis=1)
        print("\nFP taxa by the main source of their fragments, and by their fragments:")
        print(pd.crosstab(fp["main"].str[4:], fp["frag_bin"]).to_string())
        print("\nFP taxa whose fragments come from one source species: "
              f"{(fp['src_species'] == 1).mean():.3f}; median top-source share {fp['src_top_share'].median():.3f}")


def fn_anatomy(t):
    section("3. Where the reads of the FN taxa's own species went (against the TP taxa's)")
    cols = ["own_counted_here", "own_lowmapq_here", "own_counted_elsewhere", "own_lowmapq_elsewhere", "own_unaligned"]
    rows = []
    for (cls, scenario), g in t[t["class"].isin(["FN", "TP"])].groupby(["class", "scenario"]):
        n = g["own_frags"].sum()
        rows.append(dict({"class": cls, "scenario": scenario, "taxa": len(g), "own fragments": n,
                          "median counted": g["frags_counted"].median(),
                          "seeded and failed per own": round(g["seeded_failed"].sum() / n, 3) if n else np.nan},
                         **{c[4:]: round(g[c].sum() / n, 3) if n else np.nan for c in cols}))
    print(pd.DataFrame(rows).to_string(index=False))
    fn = t[t["class"] == "FN"]
    if len(fn):
        print("\nFN taxa by counted fragments:")
        print(fn["frag_bin"].value_counts().sort_index().to_string())


def separation(t):
    section("4. Separation of truth by each feature (AUC; > 0.5: higher in present taxa)")
    features = SHARES + list(RATIOS)
    called, uncalled = t[t["call"] == 1], t[t["call"] == 0]
    near = t[(t["p"] - t["knob"].astype(float)).abs() <= 0.35]
    rows = []
    for f in features:
        rows.append({"feature": f,
                     "TP vs FP": auc(called["truth"], called[f]),
                     "FN vs TN": auc(uncalled["truth"], uncalled[f]),
                     "near knob": auc(near["truth"], near[f]),
                     "near, within fragment bins": binned_auc(near, f, "frag_bin"),
                     "all, within p bins": binned_auc(t, f, "p_bin"),
                     "FP mean": called.loc[called["truth"] == 0, f].mean(),
                     "TP mean": called.loc[called["truth"] == 1, f].mean(),
                     "FN mean": uncalled.loc[uncalled["truth"] == 1, f].mean(),
                     "TN mean": uncalled.loc[uncalled["truth"] == 0, f].mean()})
    out = pd.DataFrame(rows)
    out["strength"] = (out["all, within p bins"] - 0.5).abs()
    print(out.sort_values("strength", ascending=False).drop(columns="strength").round(3).to_string(index=False))
    print(f"\n(near knob: {len(near)} rows, {near['truth'].sum()} present; own_elsewhere_share needs the read's "
          "source: an oracle, for reference only)")


def best_f1(y, s):
    order = np.argsort(-s)
    y = np.asarray(y)[order]
    tp = np.cumsum(y)
    fp = np.cumsum(1 - y)
    fn = y.sum() - tp
    f1 = 2 * tp / (2 * tp + fp + fn)
    i = int(np.argmax(f1))
    return f1[i], int(fp[i]), int(fn[i])


def incremental(t, folds, seed):
    section("5. Do the read features add to the model's score? (sample-grouped CV, rows with any fragment)")
    d = t[t["frags_any"] > 0].copy()
    d["logit_p"] = np.log(d["p"].clip(1e-6, 1 - 1e-6) / (1 - d["p"].clip(1e-6, 1 - 1e-6)))
    new = [f for f in SHARES + ["lowmapq_share", "any_not_best_share", "failed_per_counted"] if f in d]
    y, groups = d["truth"].to_numpy(), (d["set"] + "/" + d["sample"]).to_numpy()
    k = min(folds, len(set(groups)))
    if k < 2 or len(set(y)) < 2:
        print("too few samples or one class only")
        return
    variants = {"logistic(p)": None, "gbm(p)": ["logit_p"], "gbm(p + read features)": ["logit_p"] + new,
                "gbm(read features only)": new}
    scores = {name: np.zeros(len(d)) for name in variants}
    for train, test in GroupKFold(n_splits=k).split(d, y, groups):
        for name, cols in variants.items():
            if cols is None:
                m = LogisticRegression().fit(d[["logit_p"]].iloc[train], y[train])
                scores[name][test] = m.predict_proba(d[["logit_p"]].iloc[test])[:, 1]
            else:
                m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, random_state=seed)
                m.fit(d[cols].iloc[train], y[train])
                scores[name][test] = m.predict_proba(d[cols].iloc[test])[:, 1]
    rows = []
    for name, s in scores.items():
        f1, fp, fn = best_f1(y, s)
        rows.append({"model": name, "log-loss": log_loss(y, s.clip(1e-6, 1 - 1e-6)), "AUC": roc_auc_score(y, s),
                     "best F1": f1, "FP": fp, "FN": fn})
    knob = d["knob"].astype(float).to_numpy()
    called = d["p"].to_numpy() >= knob
    tp = int((called & (y == 1)).sum())
    fp_, fn_ = int((called & (y == 0)).sum()), int((~called & (y == 1)).sum())
    rows.append({"model": "the build's calls (p >= knob)", "log-loss": np.nan, "AUC": roc_auc_score(y, d["p"]),
                 "best F1": 2 * tp / (2 * tp + fp_ + fn_), "FP": fp_, "FN": fn_})
    print(pd.DataFrame(rows).round(4).to_string(index=False))
    print(f"\n{len(d)} rows of {len(set(groups))} samples; {k} folds. F1 here is over rows with fragments only "
          "(the unseen species are not rows).")
    m = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, random_state=seed)
    cols = ["logit_p"] + new
    m.fit(d[cols], y)
    from sklearn.inspection import permutation_importance
    sample = d.sample(min(len(d), 50000), random_state=seed)
    imp = permutation_importance(m, sample[cols], sample["truth"], scoring="neg_log_loss", n_repeats=3,
                                 random_state=seed)
    order = np.argsort(-imp.importances_mean)
    print("\npermutation importance (log-loss), gbm(p + read features) fitted on all rows:")
    for i in order[:15]:
        print(f"  {cols[i]:24s} {imp.importances_mean[i]:.5f}")


def read_level(folder, t):
    section("6. Read level: alignments to another species (FP taxa's reads) against alignments to the read's own "
            "(FN taxa's reads); counted records only (primary, MAPQ >= 4)")
    r = load(folder, "records")
    r = r[(r["unmapped"] == 0) & (r["primary"] == 1) & (pd.to_numeric(r["mapq"], errors="coerce") >= 4)].copy()
    r = r[r["relation"].isin(["own", "genus", "family", "other", "host"])]
    if r.empty:
        print("no counted records")
        return
    r["wrong"] = (r["relation"] != "own").astype(int)
    r["zu_per_100"] = 100 * r["zu"] / r["aln_query"].clip(lower=1)
    r["za_min_filled"] = pd.to_numeric(r["za_min"], errors="coerce").fillna(99)
    r["clipped"] = ((r["clip_left"] >= 10) | (r["clip_right"] >= 10)).astype(int)
    r["at_edge"] = ((r["at_gene_start"] == 1) | (r["at_gene_end"] == 1)).astype(int)
    r["abs_tlen"] = r["tlen"].abs()
    # The alternatives (ZA) and failed seeds (ZF) on taxa the sample calls: the read may be theirs.
    called = {(s, str(x)) for s, x in zip(t.loc[t["call"] == 1, "sample"], t.loc[t["call"] == 1, "taxon"])}

    def on_called(sample, tags, sep):
        if not isinstance(tags, str) or tags in ("", "*"):
            return 0
        return int(any((sample, part.split(":")[0]) in called for part in tags.split(sep)))
    r["za_called"] = [on_called(s, z, ",") for s, z in zip(r["sample"], r["za"])]
    if "zf" in r:
        r["zf_called"] = [on_called(s, z, ",") for s, z in zip(r["sample"], r["zf"])]
    else:
        r["zf_called"] = np.nan
    print(pd.crosstab(r["reasons"].str.split(":").str[0], r["relation"]).to_string())
    groups = {"FP-taxon records vs own alignments": r[r["reasons"].str.contains("FP:") | (r["wrong"] == 0)],
              "all records with a known source": r}
    for title, g in groups.items():
        print(f"\n{title}: {len(g)} records, {g['wrong'].sum()} to another species")
        rows = []
        for f in RECORD_FEATURES:
            if f not in g or g[f].isna().all():
                continue
            rows.append({"feature": f, "AUC (wrong)": auc(g["wrong"], g[f]),
                         "wrong mean": g.loc[g["wrong"] == 1, f].mean(), "own mean": g.loc[g["wrong"] == 0, f].mean()})
        out = pd.DataFrame(rows)
        out["strength"] = (out["AUC (wrong)"] - 0.5).abs()
        print(out.sort_values("strength", ascending=False).drop(columns="strength").round(3).to_string(index=False))


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    t = taxa_frame(opts.dir)
    overview(t)
    fp_anatomy(t)
    fn_anatomy(t)
    separation(t)
    incremental(t, opts.folds, opts.seed)
    read_level(opts.dir, t)
    return 0


if __name__ == "__main__":
    sys.exit(main())
