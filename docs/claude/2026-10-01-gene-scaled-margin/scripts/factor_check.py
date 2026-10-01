#!/usr/bin/env python3
"""factor_check.py WORLD_DIR [DB] - the gene conservation factors protal --build estimated (DB, default db7_full:
its gene_conservation.tsv, unpacked from database.protal by protal --unpack_db into DB/unpacked) against the
simulator's true rates (gtdb/simulation/gene_rates.tsv, scaled so that the median gene has rate 1, as the
factors are). Prints the correlation, the factors' range and their median by the simulator's gene category.
"""
import csv
import math
import os
import statistics
import subprocess
import sys

B = os.path.expanduser(sys.argv[1])
db = os.path.join(B, sys.argv[2] if len(sys.argv) > 2 else "db7_full")
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "..", "..", "..", "scripts", "mini_db"))
from simulate_gtdb_release import marker_category  # noqa: E402

table = os.path.join(db, "unpacked", "gene_conservation.tsv")
if not os.path.exists(table):
    protal = os.environ.get("PROTAL", os.path.expanduser("~/fix-build/bin/protal"))
    subprocess.run([protal, "--unpack_db", "--db", db, "--unpack_dir", os.path.join(db, "unpacked")],
                   check=True, stdout=subprocess.DEVNULL)
factor = {}
with open(table) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        factor[int(row["geneid"])] = float(row["factor"])
marker_of = {}
with open(os.path.join(db, "gene2geneid.tsv")) as fh:
    for line in fh:
        marker, gid = line.rstrip("\n").split("\t")
        marker_of[int(gid)] = marker
truth = {}
rates_file = os.path.join(B, "gtdb", "simulation", "gene_rates.tsv")
if os.path.exists(rates_file):
    with open(rates_file) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            truth[row["marker"]] = (float(row["rate"]), row.get("category", ""))
names = {}
with open(os.path.join(here, "..", "..", "..", "..", "scripts", "mini_db", "markers_r226.tsv")) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 3 and not line.startswith(("#", "set\t")):
            names[f[1]] = f[2]

values = sorted(factor.values())
print(f"{len(factor)} genes, factors {values[0]:.2f}-{values[-1]:.2f} (median {statistics.median(values):.2f}); "
      f"at the bounds 0.25 / 4: {sum(v <= 0.2501 for v in values)} / {sum(v >= 3.999 for v in values)}")
if truth:
    genes = [g for g in factor if marker_of.get(g) in truth]
    median_rate = statistics.median(truth[marker_of[g]][0] for g in genes)
    t = [truth[marker_of[g]][0] / median_rate for g in genes]
    e = [factor[g] for g in genes]
    mt, me = statistics.mean(t), statistics.mean(e)
    corr = sum((a - mt) * (b - me) for a, b in zip(t, e)) / math.sqrt(sum((a - mt) ** 2 for a in t) * sum((b - me) ** 2 for b in e))
    print(f"correlation with the true rates (scaled to a median of 1): {corr:.3f} over {len(genes)} genes")
    by = {}
    for g, a, b in zip(genes, t, e):
        by.setdefault(marker_category(names.get(marker_of[g], "")), []).append((a, b))
    for c, v in sorted(by.items()):
        print(f"  {c}: {len(v)} genes, true {statistics.median(a for a, _ in v):.2f}, estimated {statistics.median(b for _, b in v):.2f}")
