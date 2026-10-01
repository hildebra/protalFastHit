#!/usr/bin/env python3
"""feature_by_class.py [BENCH_DIR] - the conservation pattern of taxa by what put reads on them.

On the 0.7.1 pipeline's paired-end training table (BENCH_DIR/V071/training/training_data.tsv; default ~/bench071) and
its profiles' gene logs (VCov per hit gene, all reads the profiler keeps), the two features as protal 25d457e computes
them (Taxon::ConservationPattern): log2 of the median depth of a taxon's hit genes with factor below 1 (the 0.7.1
database's gene_conservation.tsv, unpacked by features_exp.py) over that of its other hit genes, each + 0.001, and the
conserved share of its hit genes. Groups: present taxa beside no congener held out from the database, present taxa
beside one (meta_novel_congener), and absent taxa whose closest species in the sample is a held-out congener
(meta_novel_level species, meta_relative_rank genus: they hold that species' reads), or anything else. Writes
BENCH_DIR/results_conservation/feature_by_class.md.
"""
import collections
import glob
import math
import os
import statistics
import sys

import pandas as pd

B = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else "~/bench071")
V = os.path.join(B, "V071")
OUT = os.path.join(B, "results_conservation")


def main():
    fac = pd.read_csv(os.path.join(V, "protal_db", "unpacked", "gene_conservation.tsv"), sep="\t")
    conserved = set(fac.loc[fac["factor"] < 1, "geneid"])
    pattern = {}
    for path in glob.glob(os.path.join(V, "training", "points", "rl*", "protal", "profiles", "*.profile.genes.log")):
        sample = os.path.basename(path)[:-len(".profile.genes.log")]
        t = pd.read_csv(path, sep="\t", usecols=["TaxID", "GeneID", "VCov"])
        for taxid, g in t.groupby("TaxID"):
            slow = g.loc[g["GeneID"].isin(conserved), "VCov"]
            fast = g.loc[~g["GeneID"].isin(conserved), "VCov"]
            ratio = math.log2((slow.median() + 1e-3) / (fast.median() + 1e-3)) if len(slow) and len(fast) else 0.0
            pattern[(sample, int(taxid))] = (ratio, len(slow) / len(g), len(g))
    df = pd.read_csv(os.path.join(V, "training", "training_data.tsv"), sep="\t", low_memory=False)
    groups = collections.defaultdict(list)
    for row in df.itertuples():
        p = pattern.get((row.meta_sample, int(row.taxon)))
        if p is None:
            continue
        if row.truth:
            group = "present, a congener held out in the sample" if row.meta_novel_congener else "present, no congener held out"
        elif row.meta_novel_level == "species" and row.meta_relative_rank == "genus":
            group = "absent, congener of a held-out species in the sample"
        else:
            group = "absent, other"
        groups[group].append(p)
    lines = ["# The conservation pattern by class of taxon (feature_by_class.py)", "",
             "0.7.1 pipeline, paired-end training samples. Only taxa with 5 or more hit genes in the second part.", "",
             "| taxa | rows | ratio median (all) | ratio quartiles (all) | rows with >= 5 genes | ratio median (>= 5 genes) "
             "| conserved share median (>= 5 genes) |", "|---|---|---|---|---|---|---|"]
    for group in sorted(groups):
        vals = groups[group]
        r = sorted(v[0] for v in vals)
        many = [v for v in vals if v[2] >= 5]
        q = statistics.quantiles(r, n=4) if len(r) > 3 else [float("nan")] * 3
        lines.append(f"| {group} | {len(vals)} | {statistics.median(r):+.3f} | {q[0]:+.3f} / {q[2]:+.3f} | {len(many)} | "
                     + (f"{statistics.median(v[0] for v in many):+.3f} | {statistics.median(v[1] for v in many):.3f} |"
                        if many else "- | - |"))
    lines += ["", f"genes with factor < 1: {len(conserved)} of {len(fac)}"]
    text = "\n".join(lines) + "\n"
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "feature_by_class.md"), "w") as fh:
        fh.write(text)
    print(text)


if __name__ == "__main__":
    main()
