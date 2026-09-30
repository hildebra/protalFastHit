#!/usr/bin/env python3
"""Compare an estimated (IQ-TREE) strain tree with the simulated true tree.

usage: treecmp.py <est.treefile> <true.nwk> --ref <reference row name> [--mix NAME=MAJORTIP]
Prints one JSON object with:
  n_tips, rf, max_rf, nrf           RF on sample tips only (reference pruned)
  rf_ref, max_rf_ref, nrf_ref       RF with the reference as a tip at the true root
  quartet                           quartet distance (sample tips only)
  pat_pearson, pat_spearman         patristic-distance correlation with the truth
  pat_ratio                         sum(est patristic) / sum(true patristic)
  term_ratio, int_ratio             est/true total terminal and internal branch length
  int_bl_pearson                    correlation of internal branch lengths over shared splits
  true_split_support_mean/min, n_true_splits_found, n_true_splits_ge95, n_true_splits
  false_split_support_mean, n_false_splits, n_false_ge95
  mix_sister                        tips in the smallest clade joining the mixed sample
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import phylo_utils as pu  # noqa: E402


def smallest_clade_with(tree, name):
    """Leaves of the sister group of `name` in the unrooted est tree (smallest side)."""
    sp = pu.splits(tree)
    leaves = set(pu.leaf_names(tree))
    best = None
    for s in sp:
        for side in (s, frozenset(leaves - s)):
            if name in side and len(side) >= 2:
                if best is None or len(side) < len(best):
                    best = side
    return sorted(best - {name}) if best else []


def compare(est_path, true_path, ref, mix=None):
    est = pu.fix_parents(pu.read_tree(est_path))
    true = pu.fix_parents(pu.read_tree(true_path))
    est_leaves = set(pu.leaf_names(est))
    mix_name, mix_major = (mix.split("=") if mix else (None, None))
    samples = sorted(l for l in est_leaves if l != ref and l != mix_name)
    true_leaves = set(pu.leaf_names(true))
    samples = [s for s in samples if s in true_leaves]
    out = {"n_tips": len(samples)}
    if len(samples) < 4:
        return out
    e = pu.fix_parents(pu.unroot(pu.prune(est, samples)))
    t = pu.fix_parents(pu.unroot(pu.prune(true, samples)))
    rf, mx, _, _ = pu.rf(e, t)
    out.update(rf=rf, max_rf=mx, nrf=rf / mx if mx else float("nan"))
    # with reference: true tree + ref as zero-length tip at the root
    if ref in est_leaves:
        t2 = pu.fix_parents(pu.read_tree(true_path))
        t2 = pu.prune(t2, samples)
        t2.add(pu.Node(ref, 0.0))
        pu.fix_parents(t2)
        e2 = pu.fix_parents(pu.prune(est, samples + [ref]))
        rf2, mx2, _, _ = pu.rf(e2, t2)
        out.update(rf_ref=rf2, max_rf_ref=mx2, nrf_ref=rf2 / mx2 if mx2 else float("nan"))
    out["quartet"] = pu.quartet_distance(e, t)
    pe = pu.patristic(e)
    pt = pu.patristic(t)
    pairs = [(a, b) for i, a in enumerate(samples) for b in samples[i + 1:]]
    x = [pt[p] for p in pairs]
    y = [pe[p] for p in pairs]
    out["pat_pearson"] = pu.pearson(x, y)
    out["pat_spearman"] = pu.spearman(x, y)
    out["pat_ratio"] = sum(y) / sum(x)
    te = pu.terminal_lengths(e)
    tt = pu.terminal_lengths(t)
    out["term_ratio"] = sum(te.values()) / sum(tt.values())
    se = pu.splits(e, with_info=True)
    st = pu.splits(t, with_info=True)
    int_e = sum(v[0] for v in se.values())
    int_t = sum(v[0] for v in st.values())
    out["int_ratio"] = int_e / int_t if int_t else float("nan")
    shared = [s for s in se if s in st]
    if len(shared) >= 3:
        out["int_bl_pearson"] = pu.pearson([st[s][0] for s in shared], [se[s][0] for s in shared])
    # terminal branch correlation
    out["term_bl_pearson"] = pu.pearson([tt[s] for s in samples], [te[s] for s in samples])
    # supports: computed on the tree with the reference pruned (supports of merged root edges kept)
    true_sup = [se[s][1] for s in shared if se[s][1] is not None]
    false = [s for s in se if s not in st]
    false_sup = [se[s][1] for s in false if se[s][1] is not None]
    out["n_true_splits"] = len(st)
    out["n_true_splits_found"] = len(shared)
    out["n_true_splits_ge95"] = sum(1 for v in true_sup if v >= 95)
    out["true_split_support_mean"] = sum(true_sup) / len(true_sup) if true_sup else None
    out["true_split_support_min"] = min(true_sup) if true_sup else None
    out["n_false_splits"] = len(false)
    out["n_false_ge95"] = sum(1 for v in false_sup if v >= 95)
    out["false_split_support_mean"] = sum(false_sup) / len(false_sup) if false_sup else None
    out["false_splits"] = [sorted(s) + [f"sup={se[s][1]}"] for s in false]
    # missed true splits and their true branch length
    out["missed_true_splits"] = [sorted(s) + [f"bl={st[s][0]:.2e}"] for s in st if s not in se]
    if mix_name and mix_name in est_leaves:
        e3 = pu.fix_parents(pu.prune(est, [l for l in est_leaves if l != ref]))
        sis = smallest_clade_with(e3, mix_name)
        out["mix_sister"] = sis
        out["mix_ok"] = (sis == [mix_major])
        term = {l.name: l.length for l in e3.leaves()}
        out["mix_term_len"] = term.get(mix_name)
        out["median_term_len"] = sorted(term.values())[len(term) // 2]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("est")
    ap.add_argument("true")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--mix", default=None)
    a = ap.parse_args()
    print(json.dumps(compare(a.est, a.true, a.ref, a.mix)))


if __name__ == "__main__":
    main()
