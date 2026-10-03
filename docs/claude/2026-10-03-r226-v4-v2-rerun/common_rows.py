"""Species-held-out scores of two runs on the rows both have (v4's samples): F1 at 0.5, log loss, FP, FN, overall and
by design depth. Usage: common_rows.py A_DIR B_DIR"""
import os
import sys

import numpy as np
import pandas as pd

a_dir, b_dir = sys.argv[1:3]


def stats(y, p):
    c = p >= 0.5
    tp, fp, fn = int((c & y).sum()), int((c & ~y).sum()), int((~c & y).sum())
    q = np.clip(p, 1e-6, 1 - 1e-6)
    ll = float(-np.mean(np.where(y, np.log(q), np.log(1 - q))))
    return 2 * tp / (2 * tp + fp + fn), ll, fp, fn


for t, suffix in (("pe", ""), ("se", "_se"), ("pb", "_pb"), ("ont", "_ont")):
    cols = ["meta_sample", "meta_read_pairs", "taxon", "truth", "p_species"]
    a = pd.read_csv(os.path.join(a_dir, f"trained_model{suffix}.predictions.tsv.gz"), sep="\t", usecols=cols)
    b = pd.read_csv(os.path.join(b_dir, f"trained_model{suffix}.predictions.tsv.gz"), sep="\t", usecols=cols)
    m = a.merge(b, on=["meta_sample", "meta_read_pairs", "taxon", "truth"], suffixes=("_a", "_b"))
    y = m["truth"].to_numpy() == 1
    fa, la, pa, na = stats(y, m["p_species_a"].to_numpy())
    fb, lb, pb_, nb = stats(y, m["p_species_b"].to_numpy())
    print(f"## {t}: {len(m)} common rows of {len(a)} and {len(b)} ({m.meta_sample.nunique()} samples)")
    print(f"all       F1 {fa:.4f} -> {fb:.4f}  log loss {la:.4f} -> {lb:.4f}  FP {pa} -> {pb_}  FN {na} -> {nb}")
    for depth, g in m.groupby("meta_read_pairs"):
        yy = g["truth"].to_numpy() == 1
        fa, la, pa, na = stats(yy, g["p_species_a"].to_numpy())
        fb, lb, pb_, nb = stats(yy, g["p_species_b"].to_numpy())
        print(f"{depth:>10} F1 {fa:.4f} -> {fb:.4f}  log loss {la:.4f} -> {lb:.4f}  FP {pa} -> {pb_}  FN {na} -> {nb}  "
              f"({g.meta_sample.nunique()} samples)")
