#!/usr/bin/env python3
"""Per-sample F1 of two r226 builds' presence models, from the builds' own predictions (no refit).

For each build (--builds NAME=DIR,...) and read type: the design's test set (trained_model*.test_predictions.tsv.gz, by
the final model), the scenarios' hold-out samples (*.scenario_predictions.tsv.gz, by the final model) and their
hold-in samples (*.predictions.tsv.gz, p_species: species held out), each joined row for row with its table for the
sample's depth (sample_log_fragments, log10 of the sample's fragments). Calls at the model's knob (metrics.json
global_knob, else 0.5), as protal calls; F1 at 0.5 and the best F1 at a threshold of the set's own beside it.

Writes OUT/per_sample.tsv (one row per sample and set) and OUT/per_set.tsv (each set pooled), and prints both.

    python3 per_sample.py --builds v13=local/v13,v14=local/v14 --out .
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

SUFFIX = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
TABLE = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
         "ont": "training_data_ont.tsv"}
COLS = ["meta_design", "meta_sample", "meta_scenario", "taxon", "truth", "sample_log_fragments", "fragments"]


def table(path):
    df = pd.read_csv(path, sep="\t", usecols=COLS, low_memory=False, dtype={"meta_scenario": str})
    df["meta_scenario"] = df["meta_scenario"].fillna("")
    df["truth"] = (df["truth"].astype(str).str.lower().isin(["1", "true", "present", "1.0"])).astype(int)
    return df


def knob_of(build, rt):
    with open(os.path.join(build, "model_logs", f"trained_model{SUFFIX[rt]}.metrics.json")) as fh:
        m = json.load(fh)
    k = (m.get("depth_knobs") or {}).get("global_knob")
    return float(k) if k is not None else 0.5


def attach(rows, pred, column):
    assert len(rows) == len(pred) and (rows["taxon"].to_numpy() == pred["taxon"].to_numpy()).all(), "misaligned"
    rows = rows.copy()
    rows["p"] = pred[column].to_numpy()
    return rows


def counts(y, p, t):
    call = p >= t
    tp, fp, fn = int((call & (y == 1)).sum()), int((call & (y == 0)).sum()), int((~call & (y == 1)).sum())
    return tp, fp, fn, (2 * tp / (2 * tp + fp + fn) if tp else 0.0)


def best(y, p):
    grid = np.round(np.arange(0.05, 0.96, 0.01), 2)
    f = [counts(y, p, t)[3] for t in grid]
    i = int(np.argmax(f))
    return f[i], float(grid[i])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--builds", required=True, help="NAME=DIR,... (each DIR with training/, test/, model_logs/)")
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    ap.add_argument("--out", required=True)
    opts = ap.parse_args()
    per_sample, per_set = [], []
    for spec in opts.builds.split(","):
        name, build = spec.split("=", 1)
        for rt in opts.read_types.split(","):
            knob = knob_of(build, rt)
            ml = os.path.join(build, "model_logs")
            train = table(os.path.join(build, "training", TABLE[rt]))
            test = table(os.path.join(build, "test", TABLE[rt]))
            pred = pd.read_csv(os.path.join(ml, f"trained_model{SUFFIX[rt]}.predictions.tsv.gz"), sep="\t",
                               usecols=["taxon", "p_species"])
            tpred = pd.read_csv(os.path.join(ml, f"trained_model{SUFFIX[rt]}.test_predictions.tsv.gz"), sep="\t",
                                usecols=["taxon", "p"])
            spred = pd.read_csv(os.path.join(ml, f"trained_model{SUFFIX[rt]}.scenario_predictions.tsv.gz"), sep="\t",
                                usecols=["taxon", "p"])
            hold_in = attach(train, pred, "p_species")
            hold_in = hold_in[hold_in["meta_scenario"] != ""]
            design = attach(test[test["meta_scenario"] == ""], tpred, "p")
            hold_out = attach(test[test["meta_scenario"] != ""], spred, "p")
            sets = [("design test", design.assign(meta_scenario="design"))]
            sets += [(f"{sc} hold-out", g) for sc, g in hold_out.groupby("meta_scenario")]
            sets += [(f"{sc} hold-in, species held out", g) for sc, g in hold_in.groupby("meta_scenario")]
            for label, g in sets:
                y, p = g["truth"].to_numpy(), g["p"].to_numpy()
                tp, fp, fn, f1 = counts(y, p, knob)
                bf, bt = best(y, p)
                per_set.append({"build": name, "read type": rt, "set": label, "knob": knob, "samples": g["meta_sample"].nunique(),
                                "taxa": len(g), "TP": tp, "FP": fp, "FN": fn, "F1": f1, "F1 at 0.5": counts(y, p, 0.5)[3],
                                "best F1": bf, "at": bt})
                if label == "design test":
                    continue
                for sample, s in g.groupby("meta_sample"):
                    y, p = s["truth"].to_numpy(), s["p"].to_numpy()
                    tp, fp, fn, f1 = counts(y, p, knob)
                    per_sample.append({"build": name, "read type": rt, "set": label, "sample": sample,
                                       "log10 fragments": float(s["sample_log_fragments"].iloc[0]),
                                       "present": int(y.sum()), "taxa": len(s), "TP": tp, "FP": fp, "FN": fn, "F1": f1,
                                       "F1 at 0.5": counts(y, p, 0.5)[3],
                                       "FN with 1-2 fragments": int(((p < knob) & (y == 1) & (s["fragments"] <= 2)).sum())})
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 1000)
    a, b = pd.DataFrame(per_set), pd.DataFrame(per_sample)
    a.to_csv(os.path.join(opts.out, "per_set.tsv"), sep="\t", index=False, float_format="%.5g")
    b.to_csv(os.path.join(opts.out, "per_sample.tsv"), sep="\t", index=False, float_format="%.5g")
    print(a.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print()
    print(b.sort_values(["read type", "set", "log10 fragments"]).to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
