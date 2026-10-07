#!/usr/bin/env python3
"""Why do the design's training rows prefer a high knob (pe 0.80) and its test samples 0.5? Per depth (read pairs)
present in both designs: the rows, the present taxa, what the present and absent taxa are (another genome of the
species, an in-silico strain, a novel species' congener), F1 at 0.5 and at the knob and the best threshold, for the
training rows with species held out (p_species) and the test rows scored by the final model (p).

    python3 design_vs_test.py --build local/v15 --read-type pe
"""
import argparse
import json
import os

import numpy as np
import pandas as pd

SUFFIX = {"pe": "", "se": "_se", "pb": "_pb", "ont": "_ont"}
META = ["meta_sample", "meta_read_pairs", "meta_scenario", "meta_rep_genome", "meta_insilico_strain",
        "meta_novel_level", "meta_neighbour_rank", "taxon", "truth"]


def f1(y, c):
    tp, fp, fn = int((c & y).sum()), int((c & ~y).sum()), int((~c & y).sum())
    return (2 * tp / (2 * tp + fp + fn) if tp else 0.0), fp, fn


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", required=True)
    ap.add_argument("--read-type", default="pe")
    opts = ap.parse_args()
    ml = os.path.join(opts.build, "model_logs")
    s = SUFFIX[opts.read_type]
    with open(os.path.join(ml, f"trained_model{s}.metrics.json")) as fh:
        knob = float(((json.load(fh).get("depth_knobs") or {}).get("global_knob")) or 0.5)
    tr = pd.read_csv(os.path.join(ml, f"trained_model{s}.predictions.tsv.gz"), sep="\t", usecols=META + ["p_species"],
                     dtype={"meta_scenario": str, "meta_novel_level": str, "meta_neighbour_rank": str})
    tr = tr[tr["meta_scenario"].isna()].rename(columns={"p_species": "p"})
    te = pd.read_csv(os.path.join(ml, f"trained_model{s}.test_predictions.tsv.gz"), sep="\t", usecols=META + ["p"],
                     dtype={"meta_scenario": str, "meta_novel_level": str, "meta_neighbour_rank": str})
    te = te[te["meta_scenario"].isna()]
    grid = np.round(np.arange(0.05, 0.96, 0.01), 2)
    rows = []
    for name, d in (("training, species held out", tr), ("test, final model", te)):
        for depth, g in list(d.groupby("meta_read_pairs")) + [("all", d)]:
            y, p = g["truth"].to_numpy() == 1, g["p"].to_numpy()
            pres, absn = g[y], g[~y]
            best = max((f1(y, p >= t)[0], t) for t in grid)
            rows.append({"set": name, "read pairs": depth, "samples": g["meta_sample"].nunique(), "taxa": len(g),
                         "present": int(y.sum()),
                         "present: another genome": round(float((pres["meta_rep_genome"] == 0).mean()), 3),
                         "present: in silico": round(float((pres["meta_insilico_strain"] == 1).mean()), 3),
                         "absent: novel congener": round(float(absn["meta_novel_level"].notna().mean()), 3),
                         "F1 at 0.5": f1(y, p >= 0.5)[0], "F1 at knob": f1(y, p >= knob)[0],
                         "best F1": best[0], "at": best[1]})
    pd.set_option("display.width", 250)
    pd.set_option("display.max_rows", 200)
    out = pd.DataFrame(rows)
    print(f"# {opts.read_type}, knob {knob}")
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
