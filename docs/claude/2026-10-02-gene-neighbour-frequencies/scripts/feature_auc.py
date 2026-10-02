#!/usr/bin/env python3
"""How well the gene neighbour features alone separate present from absent taxa in a build's training and test tables:
per read type, the means of present and absent taxa and the AUC (the chance that a present taxon scores higher than an
absent one; 0.5 none, below 0.5 the absent score higher), over all taxa and over those with reads across genes.

usage: feature_auc.py <build_gtdb_database.py OUTDIR> [features...]
"""

import csv
import os
import sys

TABLES = {"pe": "training_data.tsv", "se": "training_data_se.tsv", "pb": "training_data_pb.tsv",
          "ont": "training_data_ont.tsv"}
FEATURES = ["adjacent_expected_share", "adjacent_unlikely_share", "adjacent_support"]


def auc(present, absent):
    """Mann-Whitney AUC, ties counted half."""
    if not present or not absent:
        return None
    ranked = sorted([(v, 1) for v in present] + [(v, 0) for v in absent])
    rank_sum, i = 0.0, 0
    while i < len(ranked):
        j = i
        while j < len(ranked) and ranked[j][0] == ranked[i][0]:
            j += 1
        mid = (i + j + 1) / 2  # 1-based mean rank of the tie block
        rank_sum += mid * sum(1 for k in range(i, j) if ranked[k][1])
        i = j
    n1, n0 = len(present), len(absent)
    return (rank_sum - n1 * (n1 + 1) / 2) / (n1 * n0)


def main():
    out = sys.argv[1]
    features = sys.argv[2:] or FEATURES
    print("| table | read type | taxa present / absent | with links | feature | mean present / absent | AUC | AUC with links |")
    print("|---|---|---|---|---|---|---|---|")
    for part in ("training", "test"):
        for t, name in TABLES.items():
            path = os.path.join(out, part, name)
            if not os.path.isfile(path):
                continue
            with open(path) as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            truth = [r["truth"] == "1" for r in rows]
            # A taxon with reads across genes: an adjacent_support away from its prior (0.5) or an expected/unlikely share.
            linked = [float(r.get("adjacent_expected_share", 0)) + float(r.get("adjacent_unlikely_share", 0)) > 0 or
                      abs(float(r.get("adjacent_support", 0.5)) - 0.5) > 1e-12 for r in rows]
            for feature in features:
                if feature not in rows[0]:
                    continue
                values = [float(r[feature]) for r in rows]
                p = [v for v, y in zip(values, truth) if y]
                a = [v for v, y in zip(values, truth) if not y]
                pl = [v for v, y, l in zip(values, truth, linked) if y and l]
                al = [v for v, y, l in zip(values, truth, linked) if not y and l]
                f = lambda x: "" if x is None else f"{x:.3f}"
                mean = lambda xs: sum(xs) / len(xs) if xs else None
                print(f"| {part} | {t} | {len(p)} / {len(a)} | {len(pl)} / {len(al)} | `{feature}` | "
                      f"{f(mean(p))} / {f(mean(a))} | {f(auc(p, a))} | {f(auc(pl, al))} |")


if __name__ == "__main__":
    main()
