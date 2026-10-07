#!/usr/bin/env python3
"""How precise an "unknown species of genus G" call would be, and whether it would raise F1, on every row of the v15
build: each sample's taxa with reads, the model's score and knob, and the collector's meta columns.

Scores: the training rows (the design's samples and the scenarios' hold-in samples) with species held out (p_species,
as PREFIX.calls.tsv.gz), the test rows (the independent test set and the scenarios' hold-out samples) by the final
model; called at the model's knob (pe 0.8, the others 0.5). Truth for a genus-level call: meta_novel_congener, the
sample has a species of the row's genus that the database lacks (a novel congener); meta_novel_levels "species:N"
counts the sample's novel congeners.

A genus-level call is one per sample and genus, right when that genus has a novel congener in the sample. Rules:
  oracle    every FP call reported at genus level instead (what a perfect FP detector would give)
  p band    calls with p below q
  identity  calls with identity below t
  learned   a gradient-boosted classifier of "absent, beside a novel congener" on the rows' features, p and the genus's
            other rows in the sample, cross-validated by sample on the training rows, fitted on all of them for the
            test rows. Its score s: relabel (calls with s >= c become genus-level calls), add (rows not called with
            s >= c get one) or both.
Per rule and cut: species-level TP, FP, FN, F1 (a relabelled call is no longer a species call); the genus-level calls,
the share right (precision) and the novel congeners found; and an extended F1 that counts each sample's genera with a
novel congener and a row as positives. The same cuts on both sets: a cut that does well on the training rows is
judged on the test rows.

    OMP_NUM_THREADS=1 python3 unknown_genus.py --build local/v15 > unknown_genus.txt
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold

import fragments

PREFIX = {"pe": "trained_model", "se": "trained_model_se", "pb": "trained_model_pb", "ont": "trained_model_ont"}
KNOB = {"pe": 0.8, "se": 0.5, "pb": 0.5, "ont": 0.5}
KEEP = ["meta_sample", "meta_scenario", "meta_novel_levels", "meta_novel_level", "meta_relative_rank",
        "meta_novel_congener", "meta_lineage_genus", "taxon", "taxon_name", "truth"]
NOT_FEATURES = {"truth", "prediction", "probability", "taxon"}


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--build", required=True, help="the build's local folder (predictions, model_logs, tables)")
    p.add_argument("--types", default="pe,se,pb,ont")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


def scored_rows(build, rt):
    """Every row with its set, score p, the meta columns and the model's features (the tables'), samples as
    <set>:<sample>."""
    pre = PREFIX[rt]
    train = pd.read_csv(os.path.join(build, pre + ".predictions.tsv.gz"), sep="\t", low_memory=False)
    train["p"], train["set"] = train["p_species"], "training"
    tests = [pd.read_csv(os.path.join(build, pre + ".test_predictions.tsv.gz"), sep="\t", low_memory=False)]
    sc = os.path.join(build, "model_logs", pre + ".scenario_predictions.tsv.gz")
    if os.path.isfile(sc):
        tests.append(pd.read_csv(sc, sep="\t", low_memory=False))
    test = pd.concat(tests, ignore_index=True)
    test["set"] = "test"
    rows = pd.concat([train[KEEP + ["p", "set"]], test[KEEP + ["p", "set"]]], ignore_index=True)
    rows = rows[rows["p"].notna()].copy()
    rows["meta_sample"] = rows["set"] + ":" + rows["meta_sample"].astype(str)
    rows = rows.drop_duplicates(["meta_sample", "taxon"])
    table = fragments.model_tables(build, rt, None)
    feats = [c for c in table.columns if not c.startswith("meta_") and c not in NOT_FEATURES and c != "taxon_name"
             and pd.api.types.is_numeric_dtype(table[c])]
    table = table[["meta_sample", "taxon"] + feats].drop_duplicates(["meta_sample", "taxon"])
    rows = rows.merge(table, on=["meta_sample", "taxon"], how="left", validate="one_to_one")
    rows["meta_scenario"] = rows["meta_scenario"].fillna("").astype(str).replace("", "design")
    rows["genus"] = rows["meta_lineage_genus"].fillna("").astype(str)
    rows["novel"] = pd.to_numeric(rows["meta_novel_congener"], errors="coerce").fillna(0).astype(int) > 0
    rows["call"] = rows["p"] >= KNOB[rt]
    rows["y"] = (rows["truth"] == 0) & rows["novel"]
    # The genus's other rows in the sample (as a run would see them)
    rows = rows.sort_values(["meta_sample", "genus", "p"], ascending=[True, True, False]).reset_index(drop=True)
    g = rows.groupby(["meta_sample", "genus"])
    rows["genus_rows"] = g["taxon"].transform("size")
    rows["genus_called"] = g["call"].transform("sum") - rows["call"]
    rank = g.cumcount()
    rows["genus_p_rank"] = rank + 1
    top1 = g["p"].transform("first")
    top2 = rows["p"].where(rank == 1).groupby([rows["meta_sample"], rows["genus"]]).transform("max")
    rows["genus_max_p_other"] = np.where(rank == 0, top2, top1)
    if "fragments" in rows:
        rows["genus_fragment_share"] = rows["fragments"] / g["fragments"].transform("sum").replace(0, np.nan)
    feats += ["p", "genus_rows", "genus_called", "genus_p_rank", "genus_max_p_other"] + \
             (["genus_fragment_share"] if "fragments" in rows else [])
    return rows.copy(), feats


def novel_congeners(rows):
    """{set: the samples' novel congeners (meta_novel_levels species:N), summed}."""
    out = {}
    for which, sub in rows.drop_duplicates("meta_sample").groupby("set"):
        n = 0
        for text in sub["meta_novel_levels"].fillna("").astype(str):
            for item in text.split(","):
                rank, _, k = item.partition(":")
                if rank == "species" and k.isdigit():
                    n += int(k)
        out[which] = n
    return out


def score(rows, ug):
    """Species-level and genus-level outcome of genus-level calls `ug` (a mask over rows) on top of the calls."""
    truth = rows["truth"].to_numpy() == 1
    ug = ug & (rows["genus"].to_numpy() != "")  # no genus, no genus-level call
    call = rows["call"].to_numpy() & ~ug
    tp, fp, fn = int((call & truth).sum()), int((call & ~truth).sum()), int((~call & truth).sum())
    grp = pd.DataFrame({"s": rows["meta_sample"].to_numpy(), "g": rows["genus"].to_numpy(), "ug": ug,
                        "novel": rows["novel"].to_numpy(), "present": truth}).groupby(["s", "g"]).agg(
        ug=("ug", "any"), novel=("novel", "any"), present=("present", "any"))
    calls = int(grp["ug"].sum())
    right = int((grp["ug"] & grp["novel"]).sum())
    positives = int(grp["novel"].sum())
    f1 = 2 * tp / max(1, 2 * tp + fp + fn)
    tpx, fpx, fnx = tp + right, fp + calls - right, fn + positives - right
    return {"TP": tp, "FP": fp, "FN": fn, "F1": round(f1, 4), "genus calls": calls,
            "precision": round(right / calls, 3) if calls else np.nan, "novel genera found": right,
            "wrong, a known species of G present": int((grp["ug"] & ~grp["novel"] & grp["present"]).sum()),
            "of": positives, "F1 extended": round(2 * tpx / max(1, 2 * tpx + fpx + fnx), 4)}


def learned_scores(rows, feats, seed):
    """Out-of-fold scores on the training rows (folds by sample), the fit on all training rows for the test rows."""
    s = np.full(len(rows), np.nan)
    tr = (rows["set"] == "training").to_numpy()
    X = rows[feats].astype(float).to_numpy()
    y = rows["y"].to_numpy().astype(int)

    def model():
        return HistGradientBoostingClassifier(max_iter=200, learning_rate=0.08, max_leaf_nodes=31,
                                              min_samples_leaf=40, random_state=seed)
    idx = np.flatnonzero(tr)
    for a, b in GroupKFold(n_splits=5).split(idx, y[idx], rows["meta_sample"].to_numpy()[idx]):
        s[idx[b]] = model().fit(X[idx[a]], y[idx[a]]).predict_proba(X[idx[b]])[:, 1]
    te = ~tr
    if te.any():
        s[te] = model().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return s


def sweep(rows, masks, label):
    """Scores of a family of rules {cut: mask}; -> DataFrame."""
    out = []
    for cut, mask in masks.items():
        out.append({"rule": label, "cut": cut, **score(rows, mask)})
    return pd.DataFrame(out)


def section(rt, rows, feats, seed):
    knob = KNOB[rt]
    n_novel = novel_congeners(rows)
    for which in ("training", "test"):
        sub = rows[rows["set"] == which]
        fp = sub[sub["call"] & (sub["truth"] == 0)]
        print(f"{which}: {sub['meta_sample'].nunique()} samples, {len(sub)} rows; FP calls {len(fp)}: beside a novel "
              f"congener {fp['novel'].mean():.3f}, closest simulated species a novel congener "
              f"{((fp['meta_novel_level'] == 'species') & (fp['meta_relative_rank'] == 'genus')).mean():.3f}; "
              f"the samples' novel congeners {n_novel.get(which, 0)}, genera with one and a row "
              f"{sub.loc[sub['novel']].groupby(['meta_sample', 'genus']).ngroups}")
    s = learned_scores(rows, feats, seed)
    context = ["genus_rows", "genus_called", "genus_p_rank", "genus_max_p_other", "genus_fragment_share"]
    s0 = learned_scores(rows, [c for c in feats if c not in context], seed)
    rows = rows.assign(s=s, s0=s0)
    # The genus's rows in the sample: any called, fragments on the uncalled ones
    g = rows.groupby(["meta_sample", "genus"])
    rows["genus_any_called"] = g["call"].transform("any")
    rows["genus_uncalled_fragments"] = rows["fragments"].where(~rows["call"], 0).groupby(
        [rows["meta_sample"], rows["genus"]]).transform("sum")
    results = []
    for which in ("training", "test"):
        sub = rows[rows["set"] == which].reset_index(drop=True)
        call, truth = sub["call"].to_numpy(), sub["truth"].to_numpy() == 1
        none = np.zeros(len(sub), dtype=bool)
        tabs = [sweep(sub, {"-": none}, "baseline"),
                sweep(sub, {"all FP": call & ~truth, "FP beside a novel congener": call & ~truth & sub["novel"].to_numpy()},
                      "oracle"),
                sweep(sub, {q: call & (sub["p"].to_numpy() < q) for q in (0.85, 0.9, 0.95, 0.98, 0.99)
                            if q > knob}, "p band"),
                sweep(sub, {t: call & (sub["identity"].to_numpy() < t) for t in (0.95, 0.96, 0.97, 0.975, 0.98)},
                      "identity"),
                sweep(sub, {c: call & (sub["s"].to_numpy() >= c) for c in (0.3, 0.5, 0.7, 0.9)}, "learned relabel"),
                sweep(sub, {c: ~call & (sub["s"].to_numpy() >= c) for c in (0.3, 0.5, 0.7, 0.9)}, "learned add"),
                sweep(sub, {c: sub["s"].to_numpy() >= c for c in (0.3, 0.5, 0.7, 0.9)}, "learned both"),
                sweep(sub, {c: call & (sub["s0"].to_numpy() >= c) for c in (0.5,)}, "relabel, no genus context"),
                sweep(sub, {c: ~call & (sub["s0"].to_numpy() >= c) for c in (0.5,)}, "add, no genus context"),
                sweep(sub, {k: ~call & ~sub["genus_any_called"].to_numpy()
                            & (sub["genus_uncalled_fragments"].to_numpy() >= k) for k in (1, 3, 10, 30)},
                      "rule: genus without a call, fragments >=")]
        t = pd.concat(tabs, ignore_index=True)
        t.insert(0, "set", which)
        results.append(t)
        by = []
        for scen, g in sub.groupby("meta_scenario"):
            g = g.reset_index(drop=True)
            c, tr = g["call"].to_numpy(), g["truth"].to_numpy() == 1
            for label, mask in (("baseline", np.zeros(len(g), dtype=bool)), ("oracle: all FP", c & ~tr),
                                ("learned relabel 0.5", c & (g["s"].to_numpy() >= 0.5)),
                                ("learned both 0.5", g["s"].to_numpy() >= 0.5),
                                ("rule: no call, >= 3 fragments", ~c & ~g["genus_any_called"].to_numpy()
                                 & (g["genus_uncalled_fragments"].to_numpy() >= 3))):
                by.append({"set": which, "scenario": scen, "rule": label, **score(g, mask)})
        results.append(pd.DataFrame(by))
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print()
    print(results[0].to_string(index=False))
    print()
    print(results[2].to_string(index=False))
    print("\nby scenario:")
    print(pd.concat([results[1], results[3]]).to_string(index=False))


def main(argv=None):
    opts = parse_args(argv)
    for rt in opts.types.split(","):
        rows, feats = scored_rows(opts.build, rt)
        print(f"\n# {rt} (knob {KNOB[rt]}; {len(feats)} features)\n", flush=True)
        section(rt, rows, feats, opts.seed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
