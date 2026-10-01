#!/usr/bin/env python3
"""split_diag.py WORLD_DIR PROF_RULE [DB] - the premise of the gsplit rule: does a relative the database lacks make
a species' conserved genes deeper than its fast ones? For each species present in each sample (against db7_<DB>,
default missing), the median depth over its genes of all their reads (VCov of <sample>.profile.genes.log in
prof7/PROF_RULE) on genes with factor < 1 over that on genes with factor >= 1 (gene_conservation.tsv of the
database, unpacked by factor_check.py or protal --unpack_db into db7_<DB>/unpacked). Prints the quartiles of that
ratio by group: species outnumbered by a held-out congener in the sample (3 times its depth or more), species
with a held-out congener that does not outnumber it, species without one, and how many exceed 1.3 and 1.5.
"""
import collections
import csv
import glob
import os
import statistics
import sys

B = os.path.expanduser(sys.argv[1])
rule = sys.argv[2]
db = sys.argv[3] if len(sys.argv) > 3 else "missing"
factor = {}
with open(os.path.join(B, f"db7_{db}", "unpacked", "gene_conservation.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        factor[row["geneid"]] = float(row["factor"])
with open(os.path.join(B, "heldout.txt")) as fh:
    heldout = {line.strip().removeprefix("s__") for line in fh if line.strip()}
design = collections.defaultdict(list)
with open(os.path.join(B, "design", "design.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        design[row["sample"]].append(row)
genus = lambda s: s.split(" ")[0]  # noqa: E731 (species names are "Genus epithet")

groups = collections.defaultdict(list)
for sample, rows in sorted(design.items()):
    depth = collections.Counter()
    for r in rows:
        depth[r["species"]] += float(r["depth"])
    paths = glob.glob(os.path.join(B, "prof7", rule, db, sample, "*.profile.genes.log"))
    if not paths:
        continue
    genes = collections.defaultdict(lambda: ([], []))
    with open(paths[0]) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        lin, gid, vcov = header.index("Lineage"), header.index("GeneID"), header.index("VCov")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            species = f[lin].split(";s__")[-1]
            slow, fast = genes[species]
            (slow if factor.get(f[gid], 1.0) < 1 else fast).append(float(f[vcov]))
    for s in depth:
        if s in heldout and db == "missing" or s not in genes:
            continue
        slow, fast = genes[s]
        if len(slow) < 5 or len(fast) < 5 or statistics.median(fast) <= 0:
            continue
        ratio = statistics.median(slow) / statistics.median(fast)
        missing = [c for c in depth if c != s and genus(c) == genus(s) and c in heldout] if db == "missing" else []
        if any(depth[c] >= 3 * depth[s] for c in missing):
            groups["outnumbered by a held-out congener"].append(ratio)
        elif missing:
            groups["held-out congener, not outnumbered"].append(ratio)
        else:
            groups["no held-out congener"].append(ratio)
for name, values in sorted(groups.items()):
    q = statistics.quantiles(values, n=4) if len(values) > 1 else values * 3
    print(f"{name}: n={len(values)}, conserved/fast depth quartiles {q[0]:.2f} {q[1]:.2f} {q[2]:.2f}; "
          f">1.3: {sum(v > 1.3 for v in values)}, >1.5: {sum(v > 1.5 for v in values)}")
