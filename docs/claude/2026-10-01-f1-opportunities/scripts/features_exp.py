#!/usr/bin/env python3
"""features_exp.py [READ_TYPES...] - F1 of 0.7.1's presence model with candidate additions, refitted offline on the
training tables of the v0.7.1 benchmark's 0.7.1 pipeline ($V071, default ~/bench071/V071: training/ and test/).

The forest is the trainer's (64 trees, at most 128 leaves, balanced classes, the normalized features of
scripts/model_features.py at 0.7.1), as docs/claude/2026-10-01-alignment-features/exp_lib.py fits it. Scores: F1 of
species held out in 5-fold cross-validation of the training table (the trainer's 'species' scheme), and on the
independent test set (another design and seed) with a forest fitted on all training rows; seeds 1-3, the test
difference with a paired bootstrap over samples (its mean interval over seeds).

Candidates (each added to the base features alone, then the useful ones together):
  sample    the sample's context: log10 of its fragments summed over all taxa with reads, the number of those taxa,
            and the taxon's share of the fragments (the model sees a taxon's depth, not the sample's)
  congener  the taxon's genus in the sample: congeners with reads, the log ratio of the deepest congener's depth to
            the taxon's, whether the taxon is the deepest of its genus, and the congeners' share of the genus's
            fragments
  conserv   which genes its reads hit, by the genes' within-species conservation (gene_conservation.tsv of the
            database): the log ratio of the median depth of its conserved genes (factor < 1) to that of its fast
            ones, and the share of its hit genes that are conserved (a relative's reads land unevenly)
  nbr       the gene-neighbour shares protal 0.7.1 dumps but does not train on: adjacent_expected_share and
            adjacent_unlikely_share
Model variants: 256 trees; 512 leaves; scikit-learn's HistGradientBoostingClassifier. Thresholds: the knob 0.5, the
F1-optimal threshold of the cross-validated training calls, and one per depth bin (sample fragments in log10 bins)
chosen the same way. Training size: 25, 50 and 75% of the training samples.
Writes results/features_exp.md.
"""
import collections
import glob
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier

HERE = os.path.dirname(os.path.abspath(__file__))
V071 = os.path.expanduser(os.environ.get("V071", "~/bench071/V071"))
os.environ["V2"] = V071
sys.path.insert(0, os.path.join(HERE, "..", "..", "2026-10-01-alignment-features"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "..", "scripts"))
import exp_lib as X  # noqa: E402
from model_features import NORMALIZED_FEATURES as BASE  # noqa: E402

READ_TYPES = sys.argv[1:] or ["pe", "se", "pb", "ont"]
SEEDS = [int(s) for s in os.environ.get("SEEDS", "1 2 3").split()]
PROTAL = os.path.expanduser(os.environ.get("PROTAL", "~/bench071/bin/protal-0.7.1"))
OUT = os.path.join(HERE, "..", "results")


# ---- candidate features --------------------------------------------------------------------------------------

def add_sample(df):
    g = df.groupby("meta_sample")
    df["sample_log_fragments"] = np.log10(g["fragments"].transform("sum").clip(lower=1))
    df["sample_taxa"] = g["taxon"].transform("count")
    df["fragment_share"] = df["fragments"] / g["fragments"].transform("sum").clip(lower=1)
    return ["sample_log_fragments", "sample_taxa", "fragment_share"]


def add_congener(df):
    key = [df["meta_sample"], df["genus"]]
    g = df.groupby(key)
    n = g["taxon"].transform("count")
    df["congeners_with_reads"] = n - 1
    top = g["depth"].transform("max")
    # the deepest other congener: the top unless the taxon is the top, then the second
    second = g["depth"].transform(lambda s: s.nlargest(2).iloc[-1] if len(s) > 1 else 0)
    other_top = np.where(df["depth"] >= top, second, top)
    df["congener_depth_ratio"] = np.log10((other_top + 1e-3) / (df["depth"] + 1e-3))
    df["deepest_of_genus"] = (df["depth"] >= top).astype(int)
    genus_frag = g["fragments"].transform("sum").clip(lower=1)
    df["congener_fragment_share"] = 1 - df["fragments"] / genus_frag
    return ["congeners_with_reads", "congener_depth_ratio", "deepest_of_genus", "congener_fragment_share"]


def factors():
    """geneid -> conservation factor, from the 0.7.1 database (unpacked once next to it)."""
    path = os.path.join(V071, "protal_db", "unpacked", "gene_conservation.tsv")
    if not os.path.exists(path):
        subprocess.run([PROTAL, "--unpack_db", "--db", os.path.join(V071, "protal_db"), "--unpack_dir",
                        os.path.dirname(path)], check=True, stdout=subprocess.DEVNULL)
    t = pd.read_csv(path, sep="\t")
    return dict(zip(t["geneid"], t["factor"]))


def gene_logs(split, rt):
    """(sample, taxid) -> (median depth of conserved genes, of fast ones, conserved share of the hit genes)."""
    fac = factors()
    sub = "protal_se" if rt == "se" else "protal"
    out = {}
    for path in glob.glob(os.path.join(V071, split, "points", "*", sub, "profiles", "*.profile.genes.log")):
        sample = os.path.basename(path)[:-len(".profile.genes.log")]
        t = pd.read_csv(path, sep="\t", usecols=["TaxID", "GeneID", "VCov"])
        t["slow"] = t["GeneID"].map(fac).fillna(1.0) < 1
        for taxid, g in t.groupby("TaxID"):
            slow, fast = g.loc[g["slow"], "VCov"], g.loc[~g["slow"], "VCov"]
            out[(sample, int(taxid))] = (slow.median() if len(slow) else np.nan, fast.median() if len(fast) else np.nan,
                                         len(slow) / len(g))
    return out


def add_conserv(df, split, rt, cache={}):
    if (split, rt) not in cache:
        cache[(split, rt)] = gene_logs(split, rt)
    logs = cache[(split, rt)]
    vals = [logs.get((s, int(t)), (np.nan, np.nan, np.nan)) for s, t in zip(df["meta_sample"], df["taxon"])]
    slow, fast, share = (np.array(v, dtype=float) for v in zip(*vals))
    ratio = np.log2((slow + 1e-3) / (fast + 1e-3))
    df["conserved_fast_depth_ratio"] = np.where(np.isfinite(ratio), ratio, 0.0)
    df["conserved_hit_share"] = np.where(np.isfinite(share), share, 0.5)
    return ["conserved_fast_depth_ratio", "conserved_hit_share"]


def add_nbr(df):
    return ["adjacent_expected_share", "adjacent_unlikely_share"]


# ---- models and thresholds -----------------------------------------------------------------------------------

def model(kind, seed):
    if kind == "rf":
        return X.forest(seed)
    if kind == "rf256":
        return RandomForestClassifier(n_estimators=256, max_leaf_nodes=128, max_features="sqrt", class_weight="balanced",
                                      random_state=seed, n_jobs=6)
    if kind == "rf512leaves":
        return RandomForestClassifier(n_estimators=64, max_leaf_nodes=512, max_features="sqrt", class_weight="balanced",
                                      random_state=seed, n_jobs=6)
    if kind == "hgb":
        return HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, class_weight="balanced", random_state=seed)
    raise ValueError(kind)


def cv_and_test(train, test, features, kind, seed):
    from sklearn.model_selection import GroupKFold
    Xtr, ytr = train[features].to_numpy(float), train["truth"].to_numpy()
    p_cv = np.full(len(ytr), np.nan)
    for tr, te in GroupKFold(5).split(Xtr, ytr, train["taxon"].astype(str)):
        p_cv[te] = model(kind, seed).fit(Xtr[tr], ytr[tr]).predict_proba(Xtr[te])[:, 1]
    p_test = model(kind, seed).fit(Xtr, ytr).predict_proba(test[features].to_numpy(float))[:, 1]
    return p_cv, p_test


def best_threshold(y, p):
    grid = np.arange(0.05, 0.96, 0.01)
    return float(grid[np.argmax([X.scores(y, p >= t)["F1"] for t in grid])])


def depth_bins(df):
    return np.clip(np.floor(np.log10(df.groupby("meta_sample")["fragments"].transform("sum").clip(lower=1))), 2, 6).astype(int)


def per_bin_calls(train, p_cv, test, p_test):
    """Thresholds per depth bin, chosen on the cross-validated training calls; bins the training lacks keep 0.5."""
    tb, sb = depth_bins(train).to_numpy(), depth_bins(test).to_numpy()
    y = train["truth"].to_numpy()
    th = {b: best_threshold(y[tb == b], p_cv[tb == b]) for b in np.unique(tb) if (tb == b).sum() > 50}
    return np.array([p >= th.get(b, 0.5) for p, b in zip(p_test, sb)]), th


# ---- the experiment ------------------------------------------------------------------------------------------

def evaluate(rt, train, test, features, kind="rf", threshold="knob", base_calls=None, sample_frac=1.0):
    cv, tf, diffs, fps, fns = [], [], [], [], []
    for seed in SEEDS:
        tr = train
        if sample_frac < 1:
            samples = np.random.default_rng(seed).permutation(train["meta_sample"].unique())
            tr = train[train["meta_sample"].isin(samples[:max(1, int(len(samples) * sample_frac))])]
        p_cv, p_test = cv_and_test(tr, test, features, kind, seed)
        y_cv, y_test = tr["truth"].to_numpy(), test["truth"].to_numpy()
        if threshold == "knob":
            calls = p_test >= X.KNOB
        elif threshold == "cv":
            calls = p_test >= best_threshold(y_cv, p_cv)
        else:
            calls, _ = per_bin_calls(tr, p_cv, test, p_test)
        cv.append(X.scores(y_cv, p_cv >= X.KNOB)["F1"])
        s = X.scores(y_test, calls)
        tf.append(s["F1"])
        fps.append(s["FP"])
        fns.append(s["FN"])
        if base_calls is not None:
            diffs.append(X.bootstrap_diff(test, base_calls[seed], calls, n=500, seed=seed))
    d = np.array(diffs).mean(axis=0) if diffs else None
    return dict(cv=np.mean(cv), test=np.mean(tf), fp=np.mean(fps), fn=np.mean(fns), diff=d)


def main():
    lines = ["# Feature, model and threshold experiments (features_exp.py)", "",
             f"0.7.1's pipeline tables ({V071}); seeds {SEEDS}; CV = species held out (5 folds), at 0.5; test = the "
             "independent test set; delta = test F1 against the base, paired bootstrap (mean, 95% interval).", ""]
    for rt in READ_TYPES:
        train, test = X.load("training", rt), X.load("test", rt)
        extra = {}
        for name, fn in (("sample", add_sample), ("congener", add_congener), ("nbr", add_nbr)):
            extra[name] = fn(train)
            fn(test)
        extra["conserv"] = add_conserv(train, "training", rt)
        add_conserv(test, "test", rt)
        base_calls = {}
        for seed in SEEDS:
            _, p = cv_and_test(train, test, BASE, "rf", seed)
            base_calls[seed] = p >= X.KNOB
        rows = [("base (0.7.1's 28 features)", BASE, "rf", "knob", 1.0)]
        for name in ("sample", "congener", "conserv", "nbr"):
            rows.append((f"+ {name}", BASE + extra[name], "rf", "knob", 1.0))
        rows += [("base, CV-optimal threshold", BASE, "rf", "cv", 1.0),
                 ("base, threshold per depth bin", BASE, "rf", "bins", 1.0),
                 ("base, 256 trees", BASE, "rf256", "knob", 1.0),
                 ("base, 512 leaves", BASE, "rf512leaves", "knob", 1.0),
                 ("base, gradient boosting", BASE, "hgb", "knob", 1.0),
                 ("base, 25% of the training samples", BASE, "rf", "knob", 0.25),
                 ("base, 50% of the training samples", BASE, "rf", "knob", 0.5),
                 ("base, 75% of the training samples", BASE, "rf", "knob", 0.75)]
        lines += [f"## {rt} (training {len(train)} taxa, test {len(test)})", "",
                  "| variant | CV F1 | test F1 | test FP | test FN | delta (interval) |", "|---|---|---|---|---|---|"]
        for label, feats, kind, thr, frac in rows:
            r = evaluate(rt, train, test, feats, kind, thr, base_calls, frac)
            d = r["diff"]
            lines.append(f"| {label} | {r['cv']:.4f} | {r['test']:.4f} | {r['fp']:.1f} | {r['fn']:.1f} | "
                         f"{d[0]:+.4f} ({d[1]:+.4f}, {d[2]:+.4f}) |")
            print(lines[-1], flush=True)
        lines.append("")
    os.makedirs(OUT, exist_ok=True)
    name = "features_exp.md" if len(READ_TYPES) == 4 else f"features_exp_{'_'.join(READ_TYPES)}.md"
    with open(os.path.join(OUT, name), "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
