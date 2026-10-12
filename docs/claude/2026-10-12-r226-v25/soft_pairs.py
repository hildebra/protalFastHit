#!/usr/bin/env python3
"""Two builds that profiled the same samples, compared taxon by taxon (sample, taxon).

For each read type: which rows both builds have, which model features differ on the shared rows, how the
`fragments` feature moved (B - A), and the presence calls of both builds' models (species held out, and the
independent test set) on the union of their rows (a row a build lacks is not called by it). Calls that flip
between the builds are split by whether the row's `fragments` changed. F1 differences come with a paired
bootstrap over samples (both builds' calls on the same resampled samples).

Usage: soft_pairs.py BUILD_A_OUTDIR BUILD_B_OUTDIR [--read-types pe,se,pb,ont] [--boot 1000]
"""
import argparse
import gzip
import os
import re

import numpy as np
import pandas as pd

KEY = ["meta_sample", "taxon"]
KNOBS = np.round(np.arange(0.30, 0.901, 0.01), 2)


def suffix(rt):
    return "" if rt == "pe" else f"_{rt}"


def model_features(outdir, rt):
    with open(os.path.join(outdir, "model_logs", f"trained_model{suffix(rt)}.report.txt")) as handle:
        for line in handle:
            m = re.match(r"\d+ features \([^)]*\): (.*)$", line.strip())
            if m:
                return [f.strip() for f in m.group(1).split(",")]
    raise SystemExit(f"no feature list in {outdir}'s {rt} report")


def build_knob(outdir, rt):
    """The knob protal calls at: the build_metadata line's 'knob X for every sample', else 0.5."""
    path = os.path.join(outdir, "protal_db", "build_metadata.tsv")
    with open(path) as handle:
        for line in handle:
            name, _, value = line.rstrip("\n").partition("\t")
            if name == f"model_{rt}_depth_knobs":
                m = re.search(r"knob ([\d.]+) for every sample", value)
                return float(m.group(1)) if m else 0.5
    return 0.5


def load_table(outdir, rt, columns):
    path = os.path.join(outdir, "work", "training", f"training_data{suffix(rt)}.tsv")
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    use = [c for c in KEY + ["truth", "meta_rep_genome", "meta_insilico_strain"] + columns if c in header]
    return pd.read_csv(path, sep="\t", usecols=sorted(set(use)), low_memory=False)


def load_predictions(outdir, rt, kind):
    name = {"held": "predictions", "test": "test_predictions"}[kind]
    path = os.path.join(outdir, "model_logs", f"trained_model{suffix(rt)}.{name}.tsv.gz")
    col = "p_species" if kind == "held" else "p"
    frame = pd.read_csv(path, sep="\t", usecols=KEY + ["truth", col], low_memory=False)
    return frame.rename(columns={col: "p"})


def union_calls(pa, pb):
    """Rows of either build; p = -1 where a build has no row (never called)."""
    merged = pa.merge(pb, on=KEY, how="outer", suffixes=("_a", "_b"))
    merged["truth"] = merged["truth_a"].fillna(merged["truth_b"]).astype(int)
    merged["p_a"] = merged["p_a"].fillna(-1.0)
    merged["p_b"] = merged["p_b"].fillna(-1.0)
    return merged


def f1(truth, called):
    tp = np.sum(truth & called)
    fp = np.sum(~truth & called)
    fn = np.sum(truth & ~called)
    return (2 * tp / (2 * tp + fp + fn) if tp else 0.0), int(tp), int(fp), int(fn)


def best_knob(truth, p):
    scores = [f1(truth, p >= k)[0] for k in KNOBS]
    i = int(np.argmax(scores))
    return float(KNOBS[i]), scores[i]


def paired_boot(merged, knob_a, knob_b, boot, rng):
    """F1(B) - F1(A) at the given knobs, resampling samples; returns (difference, 2.5%, 97.5%, share B better)."""
    truth = merged["truth"].to_numpy().astype(bool)
    ca = merged["p_a"].to_numpy() >= knob_a
    cb = merged["p_b"].to_numpy() >= knob_b
    codes, samples = pd.factorize(merged["meta_sample"])
    n = len(samples)

    def counts(called):
        tp = np.bincount(codes, weights=(truth & called), minlength=n)
        fp = np.bincount(codes, weights=(~truth & called), minlength=n)
        fn = np.bincount(codes, weights=(truth & ~called), minlength=n)
        return tp, fp, fn

    a, b = counts(ca), counts(cb)

    def score(c, w):
        tp, fp, fn = (np.dot(w, x) for x in c)
        return 2 * tp / (2 * tp + fp + fn)

    full = score(b, np.ones(n)) - score(a, np.ones(n))
    diffs = np.empty(boot)
    for i in range(boot):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        diffs[i] = score(b, w) - score(a, w)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return full, lo, hi, float(np.mean(diffs > 0))


def flips(merged, knob_a, knob_b, changed):
    truth = merged["truth"].to_numpy().astype(bool)
    ca = merged["p_a"].to_numpy() >= knob_a
    cb = merged["p_b"].to_numpy() >= knob_b
    rows = []
    for label, mask in (("fragments changed", changed), ("fragments the same", ~changed)):
        rows.append((label, int(mask.sum()),
                     int(np.sum(mask & truth & ~ca & cb)), int(np.sum(mask & truth & ca & ~cb)),
                     int(np.sum(mask & ~truth & ca & ~cb)), int(np.sum(mask & ~truth & ~ca & cb))))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("a")
    parser.add_argument("b")
    parser.add_argument("--read-types", default="pe,se,pb,ont")
    parser.add_argument("--boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    print(f"A = {args.a}\nB = {args.b}")

    for rt in args.read_types.split(","):
        features = model_features(args.b, rt)
        if model_features(args.a, rt) != features:
            print(f"warning: {rt}: the builds' feature lists differ")
        ta = load_table(args.a, rt, features + ["fragments_all"])
        tb = load_table(args.b, rt, features + ["fragments_all"])
        shared = ta.merge(tb, on=KEY, suffixes=("_a", "_b"))
        print(f"\n# {rt}\nrows: A {len(ta)}, B {len(tb)}, both {len(shared)}; "
              f"only A {len(ta) - len(shared)} ({int(ta.merge(tb[KEY], on=KEY, how='left', indicator=True).query('_merge == \"left_only\"')['truth'].sum())} present), "
              f"only B {len(tb) - len(shared)} ({int(tb.merge(ta[KEY], on=KEY, how='left', indicator=True).query('_merge == \"left_only\"')['truth'].sum())} present)")
        truth_mismatch = int((shared["truth_a"] != shared["truth_b"]).sum())
        if truth_mismatch:
            print(f"warning: truth differs on {truth_mismatch} shared rows")

        diff = []
        for f in features + ["fragments_all"]:
            if f + "_a" not in shared or f + "_b" not in shared:
                continue
            x, y = shared[f + "_a"].to_numpy(float), shared[f + "_b"].to_numpy(float)
            differs = ~np.isclose(x, y, rtol=1e-6, atol=1e-9, equal_nan=True)
            if differs.any():
                diff.append((f, differs.mean()))
        diff.sort(key=lambda t: -t[1])
        print("features that differ on shared rows (share of rows): " +
              (", ".join(f"{f} {s:.2%}" for f, s in diff[:15]) if diff else "none") +
              (f", and {len(diff) - 15} more" if len(diff) > 15 else ""))

        fa, fb = shared["fragments_a"].to_numpy(float), shared["fragments_b"].to_numpy(float)
        allb = shared["fragments_all_b"].to_numpy(float)
        tr = shared["truth_a"].to_numpy().astype(bool)
        ch = ~np.isclose(fa, fb)
        print("fragments B - A on shared rows:")
        for label, m in (("present", tr), ("absent", ~tr)):
            d = fb[m] - fa[m]
            moved = ch[m]
            print(f"  {label:<8} rows {m.sum():>7}, changed {moved.mean():6.2%}; of the changed: median +{np.median(d[moved]) if moved.any() else 0:.2f}, "
                  f"mean +{d[moved].mean() if moved.any() else 0:.2f}; A <= 2 and B >= 3: {int(np.sum((fa[m] <= 2) & (fb[m] >= 3)))}; "
                  f"A <= 2 and fragments_all >= 3: {int(np.sum((fa[m] <= 2) & (allb[m] >= 3)))}")

        for kind in ("held", "test"):
            pa_, pb_ = load_predictions(args.a, rt, kind), load_predictions(args.b, rt, kind)
            merged = union_calls(pa_, pb_)
            merged = merged.merge(shared[KEY + ["fragments_a", "fragments_b"]], on=KEY, how="left") if kind == "held" else merged
            if kind == "held":
                changed = ~np.isclose(merged["fragments_a"].fillna(-1).to_numpy(float),
                                      merged["fragments_b"].fillna(-1).to_numpy(float))
            truth = merged["truth"].to_numpy().astype(bool)
            title = "species held out" if kind == "held" else "independent test set"
            print(f"\n{rt}, {title}: {len(merged)} rows ({truth.sum()} present) in {merged['meta_sample'].nunique()} samples")
            ka, kb = build_knob(args.a, rt), build_knob(args.b, rt)
            for label, knob_a, knob_b in (("knob 0.5", 0.5, 0.5), (f"the builds' knobs {ka}/{kb}", ka, kb)):
                sa = f1(truth, merged["p_a"].to_numpy() >= knob_a)
                sb = f1(truth, merged["p_b"].to_numpy() >= knob_b)
                d, lo, hi, share = paired_boot(merged, knob_a, knob_b, args.boot, rng)
                print(f"  {label}: A F1 {sa[0]:.4f} (TP {sa[1]}, FP {sa[2]}, FN {sa[3]}), B F1 {sb[0]:.4f} "
                      f"(TP {sb[1]}, FP {sb[2]}, FN {sb[3]}); B - A {d:+.4f} [{lo:+.4f}, {hi:+.4f}], B better in {share:.0%}")
            if kind == "held":
                bka, fka = best_knob(truth, merged["p_a"].to_numpy())
                bkb, fkb = best_knob(truth, merged["p_b"].to_numpy())
                d, lo, hi, share = paired_boot(merged, bka, bkb, args.boot, rng)
                print(f"  best knobs: A {bka} F1 {fka:.4f}, B {bkb} F1 {fkb:.4f}; B - A {d:+.4f} [{lo:+.4f}, {hi:+.4f}], "
                      f"B better in {share:.0%}")
                print("  calls that flip, A -> B at knob 0.5:   rows  FN->TP  TP->FN  FP->TN  TN->FP")
                for row in flips(merged, 0.5, 0.5, changed):
                    print(f"    {row[0]:<20} {row[1]:>9} {row[2]:>7} {row[3]:>7} {row[4]:>7} {row[5]:>7}")
                # The rows the soft count was made for: present, at most 2 fragments after the filters, 3 or more before.
                target = merged.merge(shared[KEY + ["fragments_all_b"]], on=KEY, how="left")
                m = (truth & (target["fragments_a"].fillna(99).to_numpy() <= 2)
                     & (target["fragments_all_b"].fillna(0).to_numpy() >= 3))
                ca = merged["p_a"].to_numpy() >= 0.5
                cb = merged["p_b"].to_numpy() >= 0.5
                print(f"  present, A fragments <= 2, fragments_all >= 3: {int(m.sum())} rows; FN A {int(np.sum(m & ~ca))}, "
                      f"B {int(np.sum(m & ~cb))} (knob 0.5)")
                m2 = (~truth & (target["fragments_a"].fillna(99).to_numpy() <= 2)
                      & (target["fragments_all_b"].fillna(0).to_numpy() >= 3))
                print(f"  absent, the same fragments: {int(m2.sum())} rows; FP A {int(np.sum(m2 & ca))}, B {int(np.sum(m2 & cb))}")


if __name__ == "__main__":
    main()
