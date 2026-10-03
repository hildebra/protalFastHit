"""Where a model's remaining test-set errors are, and how much F1 each kind costs. Calls at the model's knob curve
(the knob interpolated over log10 of the sample's fragments, as protal does), on the independent test set.
Also the F1 with oracle thresholds (the best one per test depth): how much better thresholds could gain.
Usage: error_budget.py RUN_DIR TEST_DIR [SUFFIX], e.g. local/v5 local/v5/test _se"""
import json
import os
import sys

import numpy as np
import pandas as pd

run, test_dir = sys.argv[1], sys.argv[2]
suffix = sys.argv[3] if len(sys.argv) > 3 else ""
metrics = json.load(open(os.path.join(run, f"trained_model{suffix}.metrics.json")))
pred = pd.read_csv(os.path.join(run, f"trained_model{suffix}.test_predictions.tsv.gz"), sep="\t", low_memory=False)
table = pd.read_csv(os.path.join(test_dir, f"training_data{suffix}.tsv"), sep="\t", low_memory=False,
                    usecols=["meta_sample", "taxon", "fragments", "identity"])
d = pred.merge(table, on=["meta_sample", "taxon"], how="left", validate="one_to_one")
curve = np.array(metrics["depth_knobs"]["curve"], dtype=float)
sample_fragments = d.groupby("meta_sample")["fragments"].transform("sum")
d["knob"] = np.interp(np.log10(np.maximum(sample_fragments, 1)), curve[:, 0], curve[:, 1])
d["call"] = d["p"] >= d["knob"]
y = d["truth"] == 1


def f1(tp, fp, fn):
    return 2 * tp / max(1, 2 * tp + fp + fn)


tp, fp, fn = int((d.call & y).sum()), int((d.call & ~y).sum()), int((~d.call & y).sum())
base = f1(tp, fp, fn)
print(f"## {run} {suffix or 'pe'}: knob curve F1 {base:.4f} (TP {tp}, FP {fp}, FN {fn}); metrics.json says "
      f"{metrics['test_depth_knobs']['F1']:.4f}")


def show(title, mask, kind):
    """Errors of `kind` (FP or FN) under `mask`, and the F1 if they were all avoided."""
    n = int(mask.sum())
    gain = f1(tp + (n if kind == "FN" else 0), fp - (n if kind == "FP" else 0), fn - (n if kind == "FN" else 0)) - base
    print(f"  {kind} {title:62} {n:5} ({n / max(1, fp if kind == 'FP' else fn):5.1%})  F1 +{gain:.4f} if avoided")


false_pos = d.call & ~y
novel = d["meta_novel_level"].fillna("").astype(str) != ""
print("false positives:")
show("closest simulated species is one the database lacks", false_pos & novel, "FP")
show("  ... and it is a congener (same genus)", false_pos & novel & (d.meta_relative_rank == "genus"), "FP")
show("others", false_pos & ~novel, "FP")
show("1 fragment", false_pos & (d.fragments <= 1), "FP")
show("2-9 fragments", false_pos & (d.fragments > 1) & (d.fragments < 10), "FP")
show(">= 10 fragments", false_pos & (d.fragments >= 10), "FP")
false_neg = ~d.call & y
strain = d["meta_rep_genome"] == 0
print("false negatives:")
show("simulated from another genome than the representative (strain)", false_neg & strain, "FN")
show("  ... with 1-10 fragments", false_neg & strain & (d.fragments <= 10), "FN")
show("the representative", false_neg & ~strain, "FN")
show("  ... with 1-10 fragments", false_neg & ~strain & (d.fragments <= 10), "FN")
show("a congener present in the sample (meta_neighbour_rank genus)", false_neg & (d.meta_neighbour_rank == "genus"),
     "FN")
show("identity below 0.97", false_neg & (d.identity < 0.97), "FN")

# Oracle thresholds: the best threshold of each test depth (an upper bound of what better thresholds could give).
grid = np.linspace(0.01, 0.99, 99)
otp = ofp = ofn = 0
for _, g in d.groupby("meta_read_pairs"):
    gy, gp = g["truth"].to_numpy() == 1, g["p"].to_numpy()
    best = max(((int((gp >= t)[gy].sum()), int((gp >= t)[~gy].sum()), int((gp < t)[gy].sum())) for t in grid),
               key=lambda c: f1(*c))
    otp, ofp, ofn = otp + best[0], ofp + best[1], ofn + best[2]
gy, gp = y.to_numpy(), d["p"].to_numpy()
glob = max(f1(int((gp >= t)[gy].sum()), int((gp >= t)[~gy].sum()), int((gp < t)[gy].sum())) for t in grid)
print(f"oracle thresholds: one for the test set F1 {glob:.4f}; one per test depth F1 {f1(otp, ofp, ofn):.4f} "
      f"(FP {ofp}, FN {ofn}); the knob curve {base:.4f}")
