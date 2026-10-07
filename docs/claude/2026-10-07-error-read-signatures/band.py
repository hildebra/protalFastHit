#!/usr/bin/env python3
"""Inside the identity band where strains and novel congeners overlap (0.95-0.985): how much of the errors live there,
what each feature separates once identity is held nearly fixed (narrow bins), and whether third_position_share (or the
other spectrum features) adds to the model's features in a boosted classifier with samples grouped. Real strains and
in-silico strains apart, since their spectra differ.

    python3 band.py --records ~/v15/err_analysis --build local/v15 > band.txt
"""
import argparse
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

import calls
import signatures

BASE = ["identity", "top_identity", "fragments", "lu_per_kb", "lu_gene_rate3", "excess_scaled_median", "low_mapq_share",
        "em_own_share", "sample_log_fragments", "cluster_genomes_log10", "mean_mapq", "hit_gene_fraction",
        "genus_spill", "relative_distance"]
SPECTRUM = ["third_position_share", "excess_conserved_fast_ratio", "gene_divergence_dispersion", "variant_sites_per_kb",
            "mate_lost_share", "failed_gene_share"]
BINS = [0.95, 0.96, 0.97, 0.98, 0.985]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--records", required=True)
    p.add_argument("--build", required=True)
    p.add_argument("--types", default="pe,se,pb,ont")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args(argv)


def cv_auc(d, cols, seed):
    X = d[cols].to_numpy(dtype=float)
    X = np.where(np.isfinite(X), X, np.nan)
    y = d["truth"].to_numpy()
    g = d["meta_sample"].to_numpy()
    s = np.zeros(len(d))
    for tr, te in GroupKFold(n_splits=5).split(X, y, g):
        m = HistGradientBoostingClassifier(max_iter=150, learning_rate=0.1, max_leaf_nodes=31, random_state=seed)
        s[te] = m.fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return roc_auc_score(y, s), average_precision_score(y, s)


def main(argv=None):
    opts = parse_args(argv)
    pd.set_option("display.width", 250)
    for rt in opts.types.split(","):
        t, features = calls.load(opts.build, rt)
        knob = calls.knob_of(opts.records, rt)
        t["class"] = calls.classes(t)
        call = t["p"] >= knob
        fp, fn = call & (t["truth"] == 0), ~call & (t["truth"] == 1)
        band = t["identity"].between(BINS[0], BINS[-1])
        print(f"\n# {rt} (knob {knob})\n")
        print(f"rows in the band {BINS[0]}-{BINS[-1]}: {int(band.sum())} of {len(t)}; of the FP {int((fp & band).sum())} "
              f"of {int(fp.sum())} ({(fp & band).sum() / max(1, fp.sum()):.2f}), of the FN {int((fn & band).sum())} of "
              f"{int(fn.sum())} ({(fn & band).sum() / max(1, fn.sum()):.2f}); band rows by class: "
              f"{t.loc[band, 'class'].value_counts().to_dict()}")
        for minf in (3, 10):
            d = t[band & (t["fragments"] >= minf) & t["class"].isin(["present strain", "present insilico", "present rep",
                                                                     "absent novel congener"])].copy()
            print(f"\n## band rows with >= {minf} fragments: {len(d)} ({d['class'].value_counts().to_dict()})\n")
            rows = []
            for col in SPECTRUM + ["identity", "lu_per_kb", "em_own_share", "excess_scaled_median", "low_mapq_share",
                                   "cluster_genomes_log10"]:
                if col not in d:
                    continue
                r = {"feature": col}
                for lo, hi in zip(BINS[:-1], BINS[1:]):
                    b = d[d["identity"].between(lo, hi)]
                    real = b[b["class"] != "present insilico"]
                    r[f"{lo}-{hi} all"] = signatures.auc(b["truth"], b[col])
                    r[f"{lo}-{hi} real strains"] = signatures.auc(real["truth"], real[col])
                rows.append(r)
            print("AUC (present high) within narrow identity bins, present strains (+ reps, + in-silico) against absent "
                  "novel congeners; 'real strains': without the in-silico strains:")
            print(pd.DataFrame(rows).round(3).to_string(index=False))
            print("  medians of third_position_share by bin and class:")
            d["bin"] = pd.cut(d["identity"], BINS)
            print(d.groupby(["bin", "class"], observed=True)["third_position_share"].agg(["size", "median"]).round(3)
                  .unstack("class").to_string())
            print("\nboosted classifier, samples grouped, 5 folds (present against absent novel congener), AUC / AP:")
            base = [c for c in BASE if c in d]
            allf = [c for c in features if c in d]
            for label, cols in (("base (identity, depth, unique k-mers, excess, MAPQ, EM, cluster, spill)", base),
                                ("base + third_position_share", base + ["third_position_share"]),
                                ("base + all spectrum features", base + [c for c in SPECTRUM if c in d]),
                                ("all model features", allf),
                                ("all model features without third_position_share",
                                 [c for c in allf if c != "third_position_share"]),
                                ("all without the spectrum features", [c for c in allf if c not in SPECTRUM])):
                a, ap = cv_auc(d, cols, opts.seed)
                print(f"  {label}: AUC {a:.4f}, AP {ap:.4f}")
            real = d[d["class"] != "present insilico"]
            a1, _ = cv_auc(real, base, opts.seed)
            a2, _ = cv_auc(real, base + ["third_position_share"], opts.seed)
            a3, _ = cv_auc(real, allf, opts.seed)
            print(f"  real strains only ({len(real)} rows): base {a1:.4f}, + third_position_share {a2:.4f}, all {a3:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
