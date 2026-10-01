#!/usr/bin/env python3
"""Post-hoc rules on the model's calls, tuned on species-held-out cross-validation of the training table and
applied unchanged to the independent test set:

  a) veto: a call whose identity is far below what present taxa have at its number of fragments
     (identity_z < -k) is dropped;
  b) genus only: an unsure call (KNOB <= p < U) is reported as its genus, not its species, when its reads fit a
     congener as well (b2: congener_share_d<delta> >= f, from the realignment with -m 3) or when a congener is a
     confident call in the same sample (b1: p >= U).

A demoted call is no species call (an absent one stops being a false positive, a present one becomes a false
negative at species level) and is scored at genus level: correct when a species of that genus was simulated in the
sample (present or held out of the database).

    python3 run_rules.py [seeds] [model: base|all]
"""
import itertools
import sys

import numpy as np
import pandas as pd

import exp_lib as L

SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 5
MODEL = sys.argv[2] if len(sys.argv) > 2 else "base"
BASE = list(L.NORMALIZED_FEATURES)
ALL = BASE + ["identity_z", "mean_mapq", "low_mapq_share", "mate_concordance", "clipped_share", "congener_share_d1",
              "other_genus_share_d1"]
pd.set_option("display.width", 230)


def tables(rt):
    train, test = L.load("training", rt), L.load("test", rt)
    model = L.identity_model(train)
    out = []
    for df, split in ((train, "training"), (test, "test")):
        df["identity_z"] = L.identity_z(df, model)
        df = df.merge(L.primary_features(split, rt), on=["meta_sample", "taxon"], how="left")
        df = df.merge(L.secondary_features(split, rt), on=["meta_sample", "taxon"], how="left")
        cols = ["mean_mapq", "low_mapq_share", "mate_concordance", "clipped_share"] + \
               [c for c in df.columns if c.startswith(("congener_share", "other_genus_share"))]
        df[cols] = df[cols].fillna(0.0)
        genera = L.simulated_genera(split)
        df["genus_simulated"] = [g in genera.get(s, set()) for s, g in zip(df["meta_sample"], df["genus"])]
        out.append(df.reset_index(drop=True))
    return out


def confident_congener(df, p, U):
    """Per row: another taxon of its genus is called with p >= U in the same sample."""
    d = pd.DataFrame({"s": df["meta_sample"], "g": df["genus"], "t": df["taxon"], "conf": p >= U})
    n = d.groupby(["s", "g"])["conf"].transform("sum")
    return (n - d["conf"].astype(int)) > 0


def apply(df, p, rule):
    """(species calls, demoted rows) of a rule."""
    call = p >= L.KNOB
    demote = np.zeros(len(df), dtype=bool)
    kind = rule["kind"]
    if kind == "veto":
        call = call & ~(df["identity_z"].to_numpy() < -rule["k"])
    elif kind in ("congener_fit", "confident_congener", "either"):
        unsure = call & (p < rule["U"])
        fit = (df[f"congener_share_d{rule['delta']}"].to_numpy() >= rule["f"]) if kind != "confident_congener" else False
        conf = confident_congener(df, p, rule["U"]).to_numpy() if kind != "congener_fit" else False
        demote = unsure & (fit | conf)
        call = call & ~demote
    return call, demote


def grid(kind):
    if kind == "veto":
        return [{"kind": kind, "k": k} for k in np.arange(0.5, 6.01, 0.25)]
    if kind == "confident_congener":
        return [{"kind": kind, "U": U} for U in (0.6, 0.7, 0.8, 0.9)]
    return [{"kind": kind, "U": U, "delta": d, "f": f} for U, d, f in
            itertools.product((0.6, 0.7, 0.8, 0.9), (0, 1, 2), np.arange(0.1, 1.01, 0.1))]


def fixed(kind):
    """The rule as the idea states it, without tuning: U = 0.8; reads 'fit as well' = within 1 edit, for most of them."""
    return {"veto": None, "confident_congener": {"kind": kind, "U": 0.8},
            "congener_fit": {"kind": kind, "U": 0.8, "delta": 1, "f": 0.5},
            "either": {"kind": kind, "U": 0.8, "delta": 1, "f": 0.5}}[kind]


def row(rt, seed, rule, params, s, y_te, demote, test, base_call, call, cv_F1):
    d = L.bootstrap_diff(test, base_call, call, n=500, seed=seed) if demote is not None else (0.0, 0.0, 0.0)
    demote = np.zeros(len(y_te), dtype=bool) if demote is None else demote
    genus_right = int(test.loc[demote, "genus_simulated"].sum())
    return {"read_type": rt, "seed": seed, "rule": rule, "params": params, "F1": s["F1"], "TP": s["TP"], "FP": s["FP"],
            "FN": s["FN"], "present": int(y_te.sum()), "demoted": int(demote.sum()),
            "demoted_absent": int((demote & (y_te == 0)).sum()), "demoted_present": int((demote & (y_te == 1)).sum()),
            "genus_right": genus_right, "dF1": d[0], "lo": d[1], "hi": d[2], "cv_F1": cv_F1}


def main():
    rows = []
    for rt in L.READ_TYPES:
        train, test = tables(rt)
        feats = ALL if MODEL == "all" else BASE
        feats = [c for c in feats if rt == "pe" or c != "mate_concordance"]
        y_tr, y_te = train["truth"].to_numpy(), test["truth"].to_numpy()
        for seed in range(1, SEEDS + 1):
            p_cv = L.cv_predict(train, feats, seed)
            p_te = L.test_predict(train, test, feats, seed)
            base_call = p_te >= L.KNOB
            rows.append(row(rt, seed, "none", "", L.scores(y_te, base_call), y_te, None, test, base_call, base_call,
                            L.scores(y_tr, p_cv >= L.KNOB)["F1"]))
            for kind in ("veto", "confident_congener", "congener_fit", "either"):
                # tuned on the training table's cross-validated calls
                best = max(grid(kind), key=lambda r: L.scores(y_tr, apply(train, p_cv, r)[0])["F1"])
                for label, rule in (("tuned", best), ("as stated", fixed(kind))):
                    if rule is None:
                        continue
                    cv_call, _ = apply(train, p_cv, rule)
                    call, demote = apply(test, p_te, rule)
                    params = ", ".join(f"{k}={v:.2f}" if isinstance(v, float) else f"{k}={v}"
                                       for k, v in rule.items() if k != "kind")
                    rows.append(row(rt, seed, f"{kind} ({label})", params, L.scores(y_te, call), y_te, demote, test,
                                    base_call, call, L.scores(y_tr, cv_call)["F1"]))
        print(rt, "done", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(f"{L.EXP}/rules_{MODEL}.tsv", sep="\t", index=False)
    # Everything reported, species and genus: how much of it is right, and how many present species are named.
    out["output_precision"] = (out["TP"] + out["genus_right"]) / (out["TP"] + out["FP"] + out["demoted"])
    out["species_sensitivity"] = out["TP"] / out["present"]
    out["genus_correct"] = np.where(out["demoted"] > 0, out["genus_right"] / out["demoted"].clip(lower=1), np.nan)
    summary = out.groupby(["read_type", "rule"], sort=False).agg(
        params=("params", lambda s: s.mode().iloc[0] if len(s) else ""), cv_F1=("cv_F1", "mean"), test_F1=("F1", "mean"),
        FP=("FP", "mean"), FN=("FN", "mean"), demoted=("demoted", "mean"), demoted_absent=("demoted_absent", "mean"),
        demoted_present=("demoted_present", "mean"), genus_correct=("genus_correct", "mean"),
        output_precision=("output_precision", "mean"), species_sensitivity=("species_sensitivity", "mean"),
        dF1=("dF1", "mean"), lo=("lo", "mean"), hi=("hi", "mean"))
    print(f"\nmodel: {MODEL} features, {SEEDS} seeds; test set F1 and its difference to no rule (bootstrap over "
          "samples); output_precision: of all reported species and genera, the share that is right\n")
    print(summary.round(4).to_string())


if __name__ == "__main__":
    main()
