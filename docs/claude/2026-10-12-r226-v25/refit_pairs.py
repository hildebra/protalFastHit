#!/usr/bin/env python3
"""Two refits of one table (refit.sh outputs) compared on the same rows, as soft_pairs.py compares two builds.

Held-out (species) and test-set F1 at knob 0.5, at each refit's chosen knob (its report) and at each one's best knob,
with a paired bootstrap over samples; the calls that flip, split by whether the row's `fragments` differs between the
two refits' training tables.

Usage: refit_pairs.py OUT_PREFIX_A OUT_PREFIX_B TABLE_A TABLE_B [--boot 1000]
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from soft_pairs import KEY, best_knob, f1, flips, paired_boot  # noqa: E402


def chosen_knob(prefix):
    """The knob the refit's report calls at ('knob X for every sample'), else 0.5."""
    with open(prefix + ".report.txt") as handle:
        for line in handle:
            m = re.search(r"knob ([\d.]+) for every sample", line)
            if m:
                return float(m.group(1))
    return 0.5


def predictions(prefix, kind):
    name, col = {"held": ("predictions", "p_species"), "test": ("test_predictions", "p")}[kind]
    frame = pd.read_csv(f"{prefix}.{name}.tsv.gz", sep="\t", usecols=KEY + ["truth", col], low_memory=False)
    return frame.rename(columns={col: "p"})


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("a")
    parser.add_argument("b")
    parser.add_argument("table_a")
    parser.add_argument("table_b")
    parser.add_argument("--boot", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    print(f"A = {args.a}\nB = {args.b}")
    ta = pd.read_csv(args.table_a, sep="\t", usecols=KEY + ["fragments"], low_memory=False)
    tb = pd.read_csv(args.table_b, sep="\t", usecols=KEY + ["fragments"], low_memory=False)
    frag = ta.merge(tb, on=KEY, suffixes=("_a", "_b"))
    ka, kb = chosen_knob(args.a), chosen_knob(args.b)
    for kind in ("held", "test"):
        merged = predictions(args.a, kind).merge(predictions(args.b, kind), on=KEY, suffixes=("_a", "_b"))
        merged["truth"] = merged["truth_a"].astype(int)
        truth = merged["truth"].to_numpy().astype(bool)
        print(f"\n{'species held out' if kind == 'held' else 'independent test set'}: {len(merged)} rows "
              f"({truth.sum()} present) in {merged['meta_sample'].nunique()} samples")
        for label, knob_a, knob_b in (("knob 0.5", 0.5, 0.5), (f"the refits' knobs {ka}/{kb}", ka, kb)):
            sa = f1(truth, merged["p_a"].to_numpy() >= knob_a)
            sb = f1(truth, merged["p_b"].to_numpy() >= knob_b)
            d, lo, hi, share = paired_boot(merged, knob_a, knob_b, args.boot, rng)
            print(f"  {label}: A F1 {sa[0]:.4f} (FP {sa[2]}, FN {sa[3]}), B F1 {sb[0]:.4f} (FP {sb[2]}, FN {sb[3]}); "
                  f"B - A {d:+.4f} [{lo:+.4f}, {hi:+.4f}], B better in {share:.0%}")
        if kind == "held":
            bka, fka = best_knob(truth, merged["p_a"].to_numpy())
            bkb, fkb = best_knob(truth, merged["p_b"].to_numpy())
            d, lo, hi, share = paired_boot(merged, bka, bkb, args.boot, rng)
            print(f"  best knobs: A {bka} F1 {fka:.4f}, B {bkb} F1 {fkb:.4f}; B - A {d:+.4f} [{lo:+.4f}, {hi:+.4f}], "
                  f"B better in {share:.0%}")
            m = merged.merge(frag, on=KEY, how="left")
            changed = ~np.isclose(m["fragments_a"].to_numpy(float), m["fragments_b"].to_numpy(float))
            print("  calls that flip, A -> B at knob 0.5:   rows  FN->TP  TP->FN  FP->TN  TN->FP")
            for row in flips(merged, 0.5, 0.5, changed):
                print(f"    {row[0]:<20} {row[1]:>9} {row[2]:>7} {row[3]:>7} {row[4]:>7} {row[5]:>7}")
            target = truth & (m["fragments_a"].to_numpy(float) <= 2) & (m["fragments_b"].to_numpy(float) >= 3)
            other = ~truth & (m["fragments_a"].to_numpy(float) <= 2) & (m["fragments_b"].to_numpy(float) >= 3)
            for label, mask, wrong in (("present, A <= 2 and B >= 3 fragments", target, "FN"),
                                       ("absent, A <= 2 and B >= 3 fragments", other, "FP")):
                for kn_label, knob_a, knob_b in (("0.5", 0.5, 0.5), ("best", bka, bkb)):
                    ca, cb = merged["p_a"].to_numpy() >= knob_a, merged["p_b"].to_numpy() >= knob_b
                    ea = np.sum(mask & (~ca if wrong == "FN" else ca))
                    eb = np.sum(mask & (~cb if wrong == "FN" else cb))
                    print(f"  {label}: {int(mask.sum())} rows; {wrong} A {int(ea)}, B {int(eb)} (knob {kn_label})")


if __name__ == "__main__":
    main()
