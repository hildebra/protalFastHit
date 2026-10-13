#!/usr/bin/env python3
"""Why a presence model rejects (or calls) a taxon: its features in a protal training dump (<profile>.truth_annotated,
written with --profile_truth), the model's probability for it (the PMML model scored as protal and the trainer do,
scripts/model_pmml.py), and what changes it: each feature in turn set to the value it has in every other taxon of a
reference group (default the sample's true positives), the probability's mean and largest change; and where the
taxon's value lies among the group's (its percentile). Needs numpy and pandas.

usage: explain_call.py MODEL.xml DUMP TAXID [--group tp|tn|all] [--top 25]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import model_pmml  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("model")
    ap.add_argument("dump")
    ap.add_argument("taxid")
    ap.add_argument("--group", choices=("tp", "tn", "all"), default="tp",
                    help="reference taxa: true positives, true negatives (then sorted by the largest fall) or all")
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()
    model = model_pmml.load_model(args.model)
    table = pd.read_csv(args.dump, sep="\t")
    row = table[table["taxon"].astype(str) == args.taxid]
    if row.empty:
        sys.exit(f"taxon {args.taxid} not in {args.dump}")
    row = row.iloc[0]
    if args.group == "tp":
        group = table[(table["truth"] == 1) & (table["prediction"] == 1)]
    elif args.group == "tn":
        group = table[(table["truth"] == 0) & (table["prediction"] == 0)]
    else:
        group = table
    group = group[group["taxon"].astype(str) != args.taxid]
    features = model.features
    x = row[features].to_frame().T.astype(float)
    p0 = float(model.predict(x)[0])
    print(f"# {args.dump}\n# taxon {args.taxid} {row['taxon_name']}: truth {row['truth']}, protal's probability "
          f"{row['probability']}, the model's here {p0:.6f}; reference group: {len(group)} {args.group} taxa")
    rows = []
    for f in features:
        values = group[f].astype(float).to_numpy()
        trial = pd.concat([x] * len(values), ignore_index=True)
        trial[f] = values
        p = model.predict(trial)
        own = float(row[f])
        below = np.sum(values < own) + 0.5 * np.sum(values == own)
        percentile = 100.0 * below / len(values) if len(values) else float("nan")
        rows.append((f, own, np.median(values) if len(values) else float("nan"), percentile, float(np.mean(p) - p0),
                     float(np.max(p) - p0)))
    rows.sort(key=lambda r: r[4] if args.group == "tn" else -r[4])
    print("feature\tvalue\tgroup_median\tpercentile_in_group\tmean_dP\tmax_dP")
    for f, own, med, pct, mean_dp, max_dp in rows[:args.top]:
        print(f"{f}\t{own:.6g}\t{med:.6g}\t{pct:.0f}\t{mean_dp:+.4f}\t{max_dp:+.4f}")
    print("# features whose value is outside the group's range (percentile 0 or 100)")
    for f, own, med, pct, mean_dp, max_dp in sorted(rows, key=lambda r: r[0]):
        if pct in (0.0, 100.0):
            print(f"{f}\t{own:.6g}\tmedian {med:.6g}\tpercentile {pct:.0f}\tmean_dP {mean_dp:+.4f}")


if __name__ == "__main__":
    main()
