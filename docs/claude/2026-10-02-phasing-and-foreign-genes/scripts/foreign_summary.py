#!/usr/bin/env python3
"""foreign_summary.py OUT ARM... - the foreign genes and calls of foreign_compare.sh's (or foreign_world.sh's) arms.

Per arm and read type, from each sample's .profile.truth_annotated (the taxa with reads: truth, call, depth) and
.profile.genes.log (LinksJudged, LinksUnlikely, Foreign per taxon and gene): the calls (TP, FP, FN, F1), the genes
judged by four or more links, the foreign ones among them in present and in absent taxa (a taxon with reads that is
not in the sample: a relative's reads, or a held-out species'), and, against the first arm, the change of the present
taxa's depth.
"""
import collections
import csv
import glob
import os
import statistics
import sys

csv.field_size_limit(1 << 30)
out, arms = sys.argv[1], sys.argv[2:]


def kind_of(sample):
    first = sample.split("_")[0]
    return "pe" if first.startswith("rl") else first


def load(arm):
    taxa, genes = {}, collections.defaultdict(list)
    for path in sorted(glob.glob(os.path.join(out, arm, "*.profile.truth_annotated"))):
        sample = os.path.basename(path).split(".profile")[0]
        for r in csv.DictReader(open(path), delimiter="\t"):
            taxa[(sample, r["taxon"])] = (int(r["truth"]), int(r["prediction"]), float(r["depth"]))
        for r in csv.DictReader(open(path.replace(".truth_annotated", ".genes.log")), delimiter="\t"):
            genes[(sample, r["TaxID"])].append((int(r["LinksJudged"]), int(r["LinksUnlikely"]), int(r["Foreign"])))
    return taxa, genes


data = {arm: load(arm) for arm in arms}
print("| arm | reads | samples | TP | FP | FN | F1 | genes judged (4+ links): present taxa | foreign | absent taxa | foreign | "
      "present taxa with a foreign gene | depth change of present taxa vs " + arms[0] + " (median, max abs) |")
print("|" + "---|" * 13)
for arm in arms:
    taxa, genes = data[arm]
    for kind in ("pe", "pb", "ont"):
        keys = [k for k in taxa if kind_of(k[0]) == kind]
        if not keys:
            continue
        samples = {k[0] for k in keys}
        tp = sum(1 for k in keys if taxa[k][0] == 1 and taxa[k][1] == 1)
        fp = sum(1 for k in keys if taxa[k][0] == 0 and taxa[k][1] == 1)
        fn = sum(1 for k in keys if taxa[k][0] == 1 and taxa[k][1] == 0)
        judged = {0: 0, 1: 0}
        foreign = {0: 0, 1: 0}
        with_foreign = 0
        for k in keys:
            truth = taxa[k][0]
            gs = genes.get(k, [])
            judged[truth] += sum(1 for j, _, _ in gs if j >= 4)
            n = sum(f for _, _, f in gs)
            foreign[truth] += n
            with_foreign += truth == 1 and n > 0
        changes = []
        base = data[arms[0]][0]
        for k in keys:
            if taxa[k][0] == 1 and k in base and base[k][2] > 0:
                changes.append(taxa[k][2] / base[k][2] - 1)
        f1 = 2 * tp / (2 * tp + fp + fn) if tp else 0
        ch = (f"{statistics.median(changes):+.4f}, {max(abs(c) for c in changes):.4f}" if changes else "-")
        print(f"| {arm} | {kind} | {len(samples)} | {tp} | {fp} | {fn} | {f1:.4f} | {judged[1]} | {foreign[1]} | "
              f"{judged[0]} | {foreign[0]} | {with_foreign} | {ch} |")
