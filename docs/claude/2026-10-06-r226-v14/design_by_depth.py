#!/usr/bin/env python3
"""The design's test set (the same rows in both builds) by depth: FP and FN of each build's final model at its knob
and at 0.5, from the builds' test predictions.

    python3 design_by_depth.py --builds v13=local/v13,v14=local/v14
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

SUFFIX = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--builds", required=True)
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    opts = ap.parse_args()
    pd.set_option("display.width", 250)
    for rt in opts.read_types.split(","):
        frames = []
        for spec in opts.builds.split(","):
            name, build = spec.split("=", 1)
            ml = os.path.join(build, "model_logs")
            with open(os.path.join(ml, f"trained_model{SUFFIX[rt]}.metrics.json")) as fh:
                knob = float(((json.load(fh).get("depth_knobs") or {}).get("global_knob")) or 0.5)
            df = pd.read_csv(os.path.join(ml, f"trained_model{SUFFIX[rt]}.test_predictions.tsv.gz"), sep="\t",
                             usecols=["meta_sample", "meta_read_pairs", "taxon", "truth", "p"])
            y, p = df["truth"].to_numpy(), df["p"].to_numpy()
            for t, label in sorted({(knob, f"{name} at {knob:g}"), (0.5, f"{name} at 0.5")}):
                c = p >= t
                frames.append(pd.DataFrame({"depth": df["meta_read_pairs"], "who": label,
                                            "FP": (c & (y == 0)).astype(int), "FN": (~c & (y == 1)).astype(int)}))
        both = pd.concat(frames)
        tab = both.groupby(["depth", "who"])[["FP", "FN"]].sum().unstack("who")
        tab.columns = [f"{b} {a}" for a, b in tab.columns]
        tab = tab[sorted(tab.columns, key=lambda c: (c.split(" ")[0], " at 0.5" in c, c.split(" ")[-1]))]
        tab.loc["all"] = tab.sum()
        print(f"== {rt} (FP and FN per depth)")
        print(tab.to_string())
        print()


if __name__ == "__main__":
    main()
