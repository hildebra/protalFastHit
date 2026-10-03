"""Headline numbers of r226 runs from their trained_model*.metrics.json: leaves, training samples and taxa; F1 at
knob 0.5 and at the knob curve, log loss, AP and false positives per sample with species held out; the same on the test
set; the learning curve's log losses (a quarter, half and all of the samples). Usage: compare_runs.py NAME=DIR ..."""
import json
import os
import sys

runs = [a.split("=", 1) for a in sys.argv[1:]]
for t, suffix in (("pe", ""), ("se", "_se"), ("pb", "_pb"), ("ont", "_ont")):
    print(f"## {t}")
    print(f"{'run':10} {'leaves':>6} {'samples':>7} {'taxa':>6} | {'F1 0.5':>7} {'F1 curve':>8} {'logloss':>7} {'AP':>6} {'FP/s':>5} | "
          f"{'test n':>6} {'F1 0.5':>7} {'F1 curve':>8} {'FP+FN curve':>11} {'logloss':>7} {'FP/s':>5} | learning curve log loss")
    for name, d in runs:
        path = os.path.join(d, f"trained_model{suffix}.metrics.json")
        if not os.path.exists(path):
            continue
        m = json.load(open(path))
        sp, te, dk, tdk = m["evaluation"]["species"], m["test"]["this one"], m["depth_knobs"], m["test_depth_knobs"]
        n_test = sum(r["samples"] for r in m["test_by_depth"])
        lc = " ".join(f"{r['log_loss']:.4f}" for r in m["learning_curve"])
        print(f"{name:10} {str(m['run']['args']['maxnodes']):>6} {m['data']['samples']:>7} {m['data']['taxa']:>6} | "
              f"{sp['F1']:.4f} {dk['F1_at_depth_knobs']:>8.4f} {sp['log_loss']:.4f} {sp['AP']:.4f} {sp['FP_per_sample']:5.2f} | "
              f"{n_test:>6} {te['F1']:.4f} {tdk['F1']:>8.4f} {int(tdk['FP']):>5}+{int(tdk['FN']):<5} {te['log_loss']:.4f} {te['FP_per_sample']:5.2f} | {lc}")
    for name, d in runs:
        path = os.path.join(d, f"trained_model{suffix}.metrics.json")
        if os.path.exists(path):
            m = json.load(open(path))
            print(f"  {name} curve: " + ",".join(f"{x}:{k}" for x, k in m["depth_knobs"]["curve"]))
