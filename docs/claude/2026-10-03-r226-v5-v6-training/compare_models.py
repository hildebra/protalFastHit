"""Headline numbers of the models of r226 runs (trained_model*.metrics.json): with species held out, F1 at knob 0.5,
AP, log loss and false positives per sample; on the independent test set, F1 and its false positives + false negatives
at knob 0.5, at the model's knob curve and at its target share of false calls (where the model has one), and log loss.
Usage: compare_models.py NAME=DIR ..."""
import json
import os
import sys

runs = [a.split("=", 1) for a in sys.argv[1:]]
for t, suffix in (("pe", ""), ("se", "_se"), ("pb", "_pb"), ("ont", "_ont")):
    print(f"## {t}")
    print(f"{'run':14} {'features':>4} | species held out: {'F1':>6} {'AP':>6} {'logloss':>7} {'FP/s':>5} | test set: "
          f"{'F1 0.5':>6} {'FP+FN':>9} | {'F1 curve':>8} {'FP+FN':>9} | {'F1 fdr':>6} {'FP+FN':>9} | {'logloss':>7} {'AP':>6}")
    for name, d in runs:
        path = os.path.join(d, f"trained_model{suffix}.metrics.json")
        if not os.path.exists(path):
            continue
        m = json.load(open(path))
        sp, te = m["evaluation"]["species"], m["test"]["this one"]
        dk, fd = m.get("test_depth_knobs"), m.get("test_false_calls")
        nf = str(m["data"].get("features", "")) or ""
        curve = f"{dk['F1']:8.4f} {int(dk['FP']):>4}+{int(dk['FN']):<4}" if dk else f"{'-':>8} {'-':>9}"
        fdr = f"{fd['F1']:6.4f} {int(fd['FP']):>4}+{int(fd['FN']):<4}" if fd else f"{'-':>6} {'-':>9}"
        print(f"{name:14} {nf:>4} | {sp['F1']:>23.4f} {sp['AP']:.4f} {sp['log_loss']:.4f} {sp['FP_per_sample']:5.2f} | "
              f"{te['F1']:>17.4f} {int(te['FP']):>4}+{int(te['FN']):<4} | {curve} | {fdr} | {te['log_loss']:.4f} {te['AP']:.4f}")
