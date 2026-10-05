#!/usr/bin/env python3
"""What the scenarios' hold-in samples did to the r226 v12 models: each read type's forest fitted with and without
them (the training table's rows with meta_scenario, or only the design's), with the scenarios' rows down-weighted, and
with a feature of the sample's complexity, the v12 trainer's settings otherwise (the feature set all four trainers
chose, 64 trees, 512 leaves for pe and se, 128 for pb and ont, balanced classes), over several seeds, and scored as
protal calls (knob 0.5; no curve, the sample's depth is a feature) on the same test tables: the design's test set (by
depth too) and each scenario's hold-out samples.

Variants (--variants):
  with scenarios            the v12 training table, as trained
  design only               its rows without meta_scenario
  scenarios x0.25           every row, the scenarios' with sample weight 0.25 (times the balanced class weights);
                            x0.1 and x0.5 the same at those weights
  with scenarios + taxa     every row, and sample_log_taxa: log10 of the taxa with reads in the sample (the sample's
                            rows), which protal knows when it profiles; not one of its features yet
  design only + taxa        the design's rows and sample_log_taxa

    python3 scenario_ablation.py --training local/v12/training --test local/v12/test --out results -t 6 --seeds 1,2,3,4,5

Writes OUT/runs.tsv (one row per read type, variant, seed and test set: TP, FP, FN, F1, rates, the best threshold),
OUT/by_depth.tsv (the design test set's errors by depth) and OUT/summary.tsv (mean and SD over seeds, and each
variant's paired difference to "with scenarios").
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import random_forest_cmdline as trainer  # noqa: E402
from model_features import DEFAULT_FEATURE_SET, feature_columns  # noqa: E402

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
LEAVES = {"pe": 512, "se": 512, "pb": 128, "ont": 128}
# name: (the scenarios' rows in, their weight, sample_log_taxa a feature)
VARIANTS = {"with scenarios": (True, 1.0, False), "design only": (False, 1.0, False),
            "scenarios x0.25": (True, 0.25, False), "with scenarios + taxa": (True, 1.0, True),
            "design only + taxa": (False, 1.0, True), "scenarios x0.1": (True, 0.1, False),
            "scenarios x0.5": (True, 0.5, False)}
BASE = "with scenarios"


def sample_log_taxa(frame):
    """log10 of each row's sample's rows: the taxa with reads in the sample."""
    return np.log10(frame.groupby("meta_sample")["truth"].transform("size").to_numpy(dtype=np.float64))


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--training", required=True)
    p.add_argument("--test", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--read-types", default="pe,se,pb,ont")
    p.add_argument("--variants", default=",".join(VARIANTS))
    p.add_argument("--seeds", default="1,2,3,4,5")
    p.add_argument("-t", "--threads", type=int, default=4)
    opts = p.parse_args()
    opts.knob = 0.5  # protal's call (trainer.scenario_row reads it)
    os.makedirs(opts.out, exist_ok=True)
    runs, depths = [], []
    for t in opts.read_types.split(","):
        began = time.time()
        train = trainer.load_table(os.path.join(opts.training, TABLES[t]))
        test = trainer.load_table(os.path.join(opts.test, TABLES[t]))
        train["sample_log_taxa"], test["sample_log_taxa"] = sample_log_taxa(train), sample_log_taxa(test)
        cols = feature_columns(train.columns, DEFAULT_FEATURE_SET)
        y = train["truth"].to_numpy()
        names = trainer.scenario_of(test)
        sets = [("design test set", names == "")] + [(f"{n} hold-out", names == n) for n in sorted(set(names) - {""})]
        design = trainer.scenario_of(train) == ""
        print(f"{t}: {len(train)} training rows, {int(design.sum())} of the design ({int(y[design].sum())} present), "
              f"{int((~design).sum())} of the scenarios ({int(y[~design].sum())} present); test sets "
              + ", ".join(f"{label} {int(mask.sum())}" for label, mask in sets), flush=True)
        for variant in opts.variants.split(","):
            scenarios_in, weight, taxa = VARIANTS[variant]
            use = cols + (["sample_log_taxa"] if taxa else [])
            X, Xt = train[use].to_numpy(dtype=np.float64), test[use].to_numpy(dtype=np.float64)
            rows = np.ones(len(y), dtype=bool) if scenarios_in else design
            weights = np.where(design, 1.0, weight)[rows]
            for seed in (int(s) for s in opts.seeds.split(",")):
                rf = RandomForestClassifier(n_estimators=64, max_leaf_nodes=LEAVES[t], min_samples_leaf=1,
                                            max_features="sqrt", class_weight="balanced", random_state=seed,
                                            n_jobs=opts.threads).fit(X[rows], y[rows], sample_weight=weights)
                prob = rf.predict_proba(Xt)[:, 1]
                for label, mask in sets:
                    frame = test[mask]
                    runs.append(trainer.scenario_row({"read type": t, "variant": variant, "seed": seed, "set": label},
                                                     frame, frame["truth"].to_numpy(), prob[mask], opts))
                mask = sets[0][1]
                call = trainer.call_scores(test[mask], prob[mask]) >= opts.knob
                truth = test.loc[mask, "truth"].to_numpy()
                for depth in sorted(set(test.loc[mask, "meta_read_pairs"])):
                    at = test.loc[mask, "meta_read_pairs"].to_numpy() == depth
                    depths.append({"read type": t, "variant": variant, "seed": seed, "depth": depth,
                                   "present": int(truth[at].sum()), "absent": int((truth[at] == 0).sum()),
                                   "FN": int(((truth == 1) & ~call & at).sum()), "FP": int(((truth == 0) & call & at).sum())})
            print(f"{t}: {variant} done, {time.time() - began:.0f} s in", flush=True)
    runs = pd.DataFrame(runs)
    runs.to_csv(os.path.join(opts.out, "runs.tsv"), sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(depths).to_csv(os.path.join(opts.out, "by_depth.tsv"), sep="\t", index=False)
    keys = ["read type", "set"]
    mean = runs.groupby(keys + ["variant"])[["F1", "FP", "FN", "best F1", "at"]].mean().unstack("variant")
    mean.columns = [f"{a} {b}" for a, b in mean.columns]
    wide = runs.pivot_table(index=keys + ["seed"], columns="variant", values="F1")
    for variant in wide.columns:
        if variant != BASE:
            d = (wide[variant] - wide[BASE]).groupby(level=keys)
            mean[f"dF1 {variant}"], mean[f"dF1 sd {variant}"] = d.mean(), d.std()
    mean.to_csv(os.path.join(opts.out, "summary.tsv"), sep="\t", float_format="%.5g")
    pd.set_option("display.width", 300)
    print(mean[[f"F1 {BASE}"] + [c for c in mean.columns if c.startswith("dF1 ") and " sd " not in c]]
          .to_string(float_format=lambda v: f"{v:+.4f}" if abs(v) < 0.5 else f"{v:.4f}"))


if __name__ == "__main__":
    main()
