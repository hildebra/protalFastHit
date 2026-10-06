#!/usr/bin/env python3
"""The soil errors split into causes (soil_errors.py --save rows), per read type, scenario and set.

False negatives, first match wins:
  1-2 fragments           too few reads to call (depth)
  divergent strain        a strain (real or in-silico) or representative whose reads differ from the reference by
                          excess_scaled_median > 0.015 (the range the congeners of novel species fill)
  beside a congener       a called present congener of the genus in the sample holds >= 5x its fragments
  other
False positives:
  near-identical relative the source (a novel species of the genus) differs by excess_scaled_median <= 0.015
  divergent relative      above it
  other                   not beside a novel species
And what F1 would be if each FN cause were removed (ceilings).

    python3 soil_partition.py --rows v13_rows.tsv.gz
"""
import argparse

import numpy as np
import pandas as pd


def genus_of(name):
    s = str(name)
    return (s[3:] if s.startswith("s__") else s).split(" ")[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rows", required=True)
    opts = ap.parse_args()
    pd.set_option("display.width", 250)
    df = pd.read_csv(opts.rows, sep="\t", low_memory=False)
    df["genus"] = df.taxon_name.map(genus_of)
    out = []
    for (rt, sc, st), g in df.groupby(["read_type", "scenario", "set"]):
        key = g.meta_sample.astype(str) + "\t" + g.genus
        top = g[g.outcome == "TP"].assign(k=key).groupby("k").fragments.max()
        congener = key.map(top).fillna(0).to_numpy()
        fn = (g.outcome == "FN").to_numpy()
        fp = (g.outcome == "FP").to_numpy()
        frag = g.fragments.to_numpy()
        ex = g.excess_scaled_median.to_numpy()
        c1 = fn & (frag < 3)
        c2 = fn & ~c1 & (ex > 0.015)
        c3 = fn & ~c1 & ~c2 & (congener >= 5 * np.maximum(frag, 1))
        c4 = fn & ~c1 & ~c2 & ~c3
        real = (g.strain == "real strain").to_numpy()
        novel = (g.meta_novel_level.astype(str) == "species").to_numpy()
        f1_ = fp & novel & (ex <= 0.015)
        f2_ = fp & novel & (ex > 0.015)
        f3_ = fp & ~novel
        tp = int((g.outcome == "TP").sum())
        nfp, nfn = int(fp.sum()), int(fn.sum())

        def f1(tp_, fp_, fn_):
            return 2 * tp_ / (2 * tp_ + fp_ + fn_)
        row = {"read type": rt, "scenario": sc, "set": st, "TP": tp, "FP": nfp, "FN": nfn, "F1": f1(tp, nfp, nfn),
               "FN 1-2 fragments": int(c1.sum()), "FN divergent": int(c2.sum()),
               "of them real strains": int((c2 & real).sum()), "FN beside congener": int(c3.sum()),
               "FN other": int(c4.sum()),
               "FP near-identical relative": int(f1_.sum()), "FP divergent relative": int(f2_.sum()),
               "FP other": int(f3_.sum()),
               "F1 without FN divergent": f1(tp + int(c2.sum()), nfp, nfn - int(c2.sum())),
               "F1 without FN 1-2 fragments": f1(tp + int(c1.sum()), nfp, nfn - int(c1.sum())),
               "F1 without FP near-identical": f1(tp, nfp - int(f1_.sum()), nfn)}
        out.append(row)
    out = pd.DataFrame(out)
    order = {"soil": 0, "soil_shallow": 1, "(design)": 2}
    out = out.sort_values(["read type", "scenario", "set"], key=lambda s: s.map(order) if s.name == "scenario" else s)
    print(out.to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
