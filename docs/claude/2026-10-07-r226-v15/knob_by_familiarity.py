#!/usr/bin/env python3
"""Why does the design's test set prefer a lower knob than the training rows? The knob is chosen on scores of
models that did not see the taxon's species (species held out); the test set is scored by the final model, which saw
most of its species in training. Test rows split by whether their species has a training row (present / any), each
part's F1 at 0.5, at the model's knob and its best threshold.

    python3 knob_by_familiarity.py --build local/v15 --read-types pe,se,pb,ont
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
    tp, fp, fn = int((c & y).sum()), int((c & ~y).sum()), int((~c & y).sum())
    return 2 * tp / (2 * tp + fp + fn) if tp else 0.0, fp, fn


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--build", required=True)
    ap.add_argument("--read-types", default="pe,se,pb,ont")
    opts = ap.parse_args()
    rows = []
    for rt in opts.read_types.split(","):
        ml = os.path.join(opts.build, "model_logs")
        with open(os.path.join(ml, f"trained_model{SUFFIX[rt]}.metrics.json")) as fh:
            knob = float(((json.load(fh).get("depth_knobs") or {}).get("global_knob")) or 0.5)
        tr = pd.read_csv(os.path.join(opts.build, "training", TABLE[rt]), sep="\t", usecols=["taxon", "truth"])
        seen_present = set(tr.loc[tr["truth"] == 1, "taxon"])
        seen_any = set(tr["taxon"])
        te = pd.read_csv(os.path.join(ml, f"trained_model{SUFFIX[rt]}.test_predictions.tsv.gz"), sep="\t",
                         usecols=["taxon", "truth", "p"])
        y, p = te["truth"].to_numpy() == 1, te["p"].to_numpy()
        groups = {"all": np.ones(len(te), bool),
                  "species present in training": te["taxon"].isin(seen_present).to_numpy(),
                  "species only absent in training": (te["taxon"].isin(seen_any) & ~te["taxon"].isin(seen_present)).to_numpy(),
                  "species not in training": (~te["taxon"].isin(seen_any)).to_numpy()}
        grid = np.round(np.arange(0.05, 0.96, 0.01), 2)
        for name, m in groups.items():
            if not m.any():
                continue
            best = max(((f1(y[m], p[m] >= t)[0], t) for t in grid))
            a, afp, afn = f1(y[m], p[m] >= 0.5)
            k, kfp, kfn = f1(y[m], p[m] >= knob)
            rows.append({"read type": rt, "test rows": name, "taxa": int(m.sum()), "present": int(y[m].sum()),
                         "F1 at 0.5": a, "FP, FN at 0.5": f"{afp}, {afn}", f"F1 at knob": k, "knob": knob,
                         "FP, FN at knob": f"{kfp}, {kfn}", "best F1": best[0], "at": best[1]})
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
