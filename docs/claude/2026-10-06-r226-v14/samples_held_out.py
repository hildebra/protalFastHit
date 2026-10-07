#!/usr/bin/env python3
"""The scenarios' hold-in samples scored by models that did not see the sample (p_samples) against models that did not
see the species (p_species), per sample with its depth: does a sample's depth still identify it (or its scenario)?

    python3 samples_held_out.py --build local/v14 --read-types pe,se,pb,ont
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

SUFFIX = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
TABLE = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
         "ont": "training_data_ont.tsv"}


def f1(y, c):
    tp, fp, fn = int((c & (y == 1)).sum()), int((c & (y == 0)).sum()), int((~c & (y == 1)).sum())
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0, fp, fn


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", required=True)
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    opts = ap.parse_args()
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 500)
    for rt in opts.read_types.split(","):
        ml = os.path.join(opts.build, "model_logs")
        with open(os.path.join(ml, f"trained_model{SUFFIX[rt]}.metrics.json")) as fh:
            knob = float(((json.load(fh).get("depth_knobs") or {}).get("global_knob")) or 0.5)
        t = pd.read_csv(os.path.join(opts.build, "training", TABLE[rt]), sep="\t",
                        usecols=["meta_sample", "meta_scenario", "taxon", "truth", "sample_log_fragments"],
                        dtype={"meta_scenario": str})
        pr = pd.read_csv(os.path.join(ml, f"trained_model{SUFFIX[rt]}.predictions.tsv.gz"), sep="\t",
                         usecols=["taxon", "p_samples", "p_species"])
        assert (t["taxon"].to_numpy() == pr["taxon"].to_numpy()).all()
        t["p_samples"], t["p_species"] = pr["p_samples"].to_numpy(), pr["p_species"].to_numpy()
        t = t[t["meta_scenario"].notna()]
        rows = []
        for (sc, s), g in t.groupby(["meta_scenario", "meta_sample"]):
            y = g["truth"].to_numpy()
            a, afp, afn = f1(y, g["p_samples"].to_numpy() >= knob)
            b, bfp, bfn = f1(y, g["p_species"].to_numpy() >= knob)
            rows.append({"scenario": sc, "sample": s.rsplit("_s_", 1)[-1], "log10 fragments": g["sample_log_fragments"].iloc[0],
                         "samples held out F1": a, "FP": afp, "FN": afn, "species held out F1": b, "FP ": bfp, "FN ": bfn,
                         "median p of present (samples / species)": f"{np.median(g['p_samples'][y == 1]):.3f} / "
                                                                    f"{np.median(g['p_species'][y == 1]):.3f}"})
        print(f"== {rt}, knob {knob}")
        print(pd.DataFrame(rows).sort_values(["scenario", "log10 fragments"]).to_string(index=False, float_format=lambda v: f"{v:.4f}"))
        print()


if __name__ == "__main__":
    main()
