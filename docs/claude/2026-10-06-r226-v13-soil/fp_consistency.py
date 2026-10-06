#!/usr/bin/env python3
"""How consistent are the false-positive taxa, and could a per-taxon prior catch them? (r226 v13 predictions)

From a build's model_logs (no tables needed): the training rows scored with species held out
(trained_model*.predictions.tsv.gz: the design's samples and the scenarios' hold-in ones) and the scenarios' hold-out
samples (*.scenario_predictions.tsv.gz). Calls at the model's knob.

1. Recurrence: per environment, how many false-positive taxa are false positives in 1, 2, 3+ samples; the false
   positive rate of a taxon's absent rows given that it was a false positive in another sample of the environment,
   against the rate otherwise.
2. Across environments: a taxon's propensity from the design's training samples only (its absent rows' false
   positive rate and mean score, with species held out; independent communities), and how well it separates the
   soil samples' false positives from their true positives among the called taxa (AUC), and what vetoing the calls
   of high-propensity taxa does to soil's F1 (the propensity cut chosen on the soil hold-in samples, applied to the
   hold-out ones). The build's held-out species are the same in both, so this is an upper bound for a prior learned
   from simulations: it knows which relatives the database lacks.

    python3 fp_consistency.py --build local/v13
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

READ_TYPES = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}


def load(build, rt):
    s = READ_TYPES[rt]
    cols = ["meta_sample", "meta_scenario", "meta_novel_level", "taxon", "taxon_name", "truth"]
    tr = pd.read_csv(os.path.join(build, "model_logs", f"trained_model{s}.predictions.tsv.gz"), sep="\t",
                     usecols=cols + ["p_species"], dtype={"meta_scenario": str, "meta_novel_level": str})
    tr = tr.rename(columns={"p_species": "p"})
    tr["set"] = "held-in"
    ho = pd.read_csv(os.path.join(build, "model_logs", f"trained_model{s}.scenario_predictions.tsv.gz"), sep="\t",
                     usecols=cols + ["p"], dtype={"meta_scenario": str, "meta_novel_level": str})
    ho["set"] = "hold-out"
    df = pd.concat([tr, ho], ignore_index=True)
    df["env"] = df.meta_scenario.fillna("").replace("", "design")
    with open(os.path.join(build, "model_logs", f"trained_model{s}.metrics.json")) as f:
        k = (json.load(f).get("depth_knobs") or {}).get("global_knob")
    df["knob"] = float(k) if k is not None else 0.5
    df["call"] = df.p >= df.knob
    df["fp"] = df.call & (df.truth == 0)
    df["tp"] = df.call & (df.truth == 1)
    df["fn"] = ~df.call & (df.truth == 1)
    return df


def f1(g, veto=None):
    call = g.call.to_numpy() if veto is None else g.call.to_numpy() & ~veto
    y = g.truth.to_numpy()
    tp, fp, fn = (call & (y == 1)).sum(), (call & (y == 0)).sum(), (~call & (y == 1)).sum()
    return 2 * tp / (2 * tp + fp + fn), int(fp), int(fn)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", required=True)
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    opts = ap.parse_args()
    pd.set_option("display.width", 250)
    for rt in opts.read_types.split(","):
        df = load(opts.build, rt)
        print(f"\n######## {rt} (knob {df.knob.iloc[0]})")
        for env in ["design", "soil", "soil_shallow", "gut"]:
            g = df[df.env == env]
            if not len(g):
                continue
            absent = g[g.truth == 0]
            per = absent.groupby("taxon").agg(rows=("fp", "size"), fp=("fp", "sum"))
            fpt = per[per.fp > 0]
            # the FP rate of a taxon's absent rows given it was FP in another sample
            a = absent.join(per.fp.rename("fp_taxon"), on="taxon")
            other = a.fp_taxon - a.fp.astype(int)
            rate_given = a.fp[other > 0].mean()
            rate_not = a.fp[other == 0].mean()
            share_rec = fpt.fp[fpt.fp >= 2].sum() / max(1, fpt.fp.sum())
            print(f"{env}: {(g.set + g.meta_sample).nunique()} samples, {int(per.fp.sum())} FP of {len(fpt)} taxa "
                  f"(1 sample {int((fpt.fp == 1).sum())}, 2: {int((fpt.fp == 2).sum())}, 3+: {int((fpt.fp >= 3).sum())}); "
                  f"FP in taxa that recur: {share_rec:.1%}; FP rate of an absent row: {absent.fp.mean():.4f}, "
                  f"if the taxon was FP in another sample: {rate_given:.4f}, else {rate_not:.4f}; "
                  f"absent rows near a novel species: {(absent.meta_novel_level.astype(str) == 'species').mean():.1%}")
        # propensity from the design's training samples only
        des = df[(df.env == "design") & (df.truth == 0)]
        prop = des.groupby("taxon").agg(n=("fp", "size"), fp=("fp", "sum"), mean_p=("p", "mean"))
        prior_rate = des.fp.mean()
        prop["rate"] = (prop.fp + 2 * prior_rate) / (prop.n + 2)  # smoothed
        for env in ["soil", "soil_shallow"]:
            g = df[df.env == env].join(prop[["rate", "mean_p", "n"]], on="taxon")
            g["rate"] = g.rate.fillna(prior_rate)
            g["mean_p"] = g.mean_p.fillna(des.p.mean())
            called = g[g.call]
            if not len(called):
                continue
            cov = called.n.notna().mean()
            y = called.fp.astype(int).to_numpy()
            auc_r = roc_auc_score(y, called.rate) if 0 < y.sum() < len(y) else float("nan")
            auc_p = roc_auc_score(y, called.mean_p) if 0 < y.sum() < len(y) else float("nan")
            hin, hout = g[g.set == "held-in"], g[g.set == "hold-out"]
            # veto calls of taxa whose design-sample FP rate is above t, t chosen on the hold-in samples
            best_t, best_f = None, f1(hin)[0]
            for t in np.quantile(prop.rate, [0.5, 0.75, 0.9, 0.95, 0.98, 0.99, 0.995]):
                f = f1(hin, (hin.rate > t).to_numpy())[0]
                if f > best_f:
                    best_t, best_f = t, f
            base = f1(hout)
            vetoed = f1(hout, (hout.rate > best_t).to_numpy()) if best_t is not None else base
            print(f"{env}: called taxa seen absent in the design's samples: {cov:.1%}; AUC of the design-sample "
                  f"propensity for FP among the called: FP rate {auc_r:.3f}, mean score {auc_p:.3f}; veto above "
                  f"{best_t if best_t is None else round(best_t, 4)} (chosen on hold-in): hold-out F1 {base[0]:.4f} "
                  f"(FP {base[1]}, FN {base[2]}) -> {vetoed[0]:.4f} (FP {vetoed[1]}, FN {vetoed[2]})")
        # soil hold-in -> soil hold-out: the same environment's other samples (shares sources: an upper bound)
        for env in ["soil", "soil_shallow"]:
            hin = df[(df.env == env) & (df.set == "held-in")]
            hout = df[(df.env == env) & (df.set == "hold-out")]
            if not len(hout):
                continue
            pin = hin[hin.truth == 0].groupby("taxon").fp.mean()
            pres_in = hin[hin.truth == 1].groupby("taxon").size()
            h = hout.join(pin.rename("fp_in"), on="taxon").join(pres_in.rename("present_in"), on="taxon")
            called = h[h.call]
            y = called.fp.astype(int).to_numpy()
            x = called.fp_in.fillna(0).to_numpy()
            auc = roc_auc_score(y, x)
            veto = (h.fp_in.fillna(0) >= 1 / 3).to_numpy()
            base, v = f1(h), f1(h, veto)
            print(f"{env} hold-in -> hold-out: called hold-out taxa FP in a hold-in sample: "
                  f"{(called.fp_in.fillna(0) > 0).mean():.1%}; AUC of the hold-in FP rate for FP among the called "
                  f"{auc:.3f}; veto if FP in >= 1/3 of its hold-in absent rows: F1 {base[0]:.4f} -> {v[0]:.4f} "
                  f"(FP {base[1]} -> {v[1]}, FN {base[2]} -> {v[2]})")


if __name__ == "__main__":
    main()
