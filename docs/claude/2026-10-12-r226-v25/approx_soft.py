#!/usr/bin/env python3
"""The soft fragment count from a table's own columns, checked where the true one is known, or written into a table.

protal's soft count (MicrobialProfile::SoftFragments) is kept + share x (all - kept), with share = (em_all x records -
em_kept x kept_records) / (records - kept_records), clamped to [0, 1]. The training tables lack the record counts; with
records in proportion to fragments it becomes share = (em_all x all - em_kept x kept) / (all - kept), from the columns
fragments (kept), fragments_all, em_own_share and em_kept_own_share.

check HARD_TABLE SOFT_TABLE: the approximation from HARD_TABLE's columns against SOFT_TABLE's fragments (the same rows,
a build with --fragment-count mapq and one with soft; long reads, whose other features do not change).
write HARD_TABLE OUT_TABLE: HARD_TABLE with fragments replaced by the approximation.
"""
import argparse

import numpy as np
import pandas as pd

KEY = ["meta_sample", "taxon"]


def approx(frame):
    kept = frame["fragments"].to_numpy(float)
    total = frame["fragments_all"].to_numpy(float)
    em_all = frame["em_own_share"].to_numpy(float)
    em_kept = frame["em_kept_own_share"].to_numpy(float)
    dropped = np.maximum(0.0, total - kept)
    with np.errstate(divide="ignore", invalid="ignore"):
        share = np.where(dropped > 0, (em_all * total - em_kept * kept) / np.where(dropped > 0, dropped, 1), em_all)
    share = np.clip(share, 0.0, 1.0)
    return np.where(dropped > 0, kept + share * dropped, kept)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("mode", choices=("check", "write"))
    parser.add_argument("hard")
    parser.add_argument("other")
    args = parser.parse_args()
    if args.mode == "check":
        cols = KEY + ["truth", "fragments", "fragments_all", "em_own_share", "em_kept_own_share"]
        hard = pd.read_csv(args.hard, sep="\t", usecols=cols, low_memory=False)
        soft = pd.read_csv(args.other, sep="\t", usecols=KEY + ["fragments"], low_memory=False)
        both = hard.merge(soft, on=KEY, suffixes=("", "_soft"))
        est = approx(both)
        true = both["fragments_soft"].to_numpy(float)
        moved = ~np.isclose(both["fragments"].to_numpy(float), true)
        err = est - true
        print(f"{len(both)} rows, {moved.sum()} with soft != kept")
        for label, m in (("all rows", np.ones(len(both), bool)), ("rows the soft count moved", moved)):
            e = err[m]
            print(f"  {label}: |error| < 0.01 in {np.mean(np.abs(e) < 0.01):.1%}, < 0.5 in {np.mean(np.abs(e) < 0.5):.1%}; "
                  f"median error {np.median(e):+.3f}, mean {e.mean():+.3f}; Spearman-like rank r "
                  f"{np.corrcoef(pd.Series(est[m]).rank(), pd.Series(true[m]).rank())[0, 1]:.4f}")
        for t, label in ((1, "present"), (0, "absent")):
            m = moved & (both["truth"].to_numpy() == t)
            print(f"  {label}, moved rows {m.sum()}: true soft - kept mean {np.mean(true[m] - both['fragments'].to_numpy(float)[m]):.2f}, "
                  f"approximation's {np.mean(est[m] - both['fragments'].to_numpy(float)[m]):.2f}")
    else:
        frame = pd.read_csv(args.hard, sep="\t", low_memory=False, dtype=str)
        numeric = frame[["fragments", "fragments_all", "em_own_share", "em_kept_own_share"]].astype(float)
        frame["fragments"] = [f"{v:.6g}" for v in approx(numeric)]
        frame.to_csv(args.other, sep="\t", index=False)
        print(f"wrote {args.other}: {len(frame)} rows")


if __name__ == "__main__":
    main()
