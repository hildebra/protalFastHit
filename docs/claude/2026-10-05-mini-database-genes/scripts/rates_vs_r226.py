#!/usr/bin/env python3
"""The mini world's per-gene evolution speeds (simulate_gtdb_release.py --gene_rates categories) against the real
GTDB r226 genes' conservation factors (gene_congeners.tsv of the r226 v10 build: within_factor, protal's factor;
between_factor, the same between congeneric species).

Gene ids follow the converter's order (bac120 markers sorted, then the ar53-only ones), which a mini database's
gene2geneid.tsv reproduces from the same markers_r226.tsv.

Usage: rates_vs_r226.py MARKERS_TSV GENE2GENEID_TSV GENE_CONGENERS_TSV
"""
import statistics
import sys

sys.path.insert(0, sys.argv[4] if len(sys.argv) > 4 else ".")
from simulate_gtdb_release import CATEGORY_RATE, marker_category  # noqa: E402


def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    r = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2
        i = j + 1
    return r


def spearman(a, b):
    ra, rb = ranks(a), ranks(b)
    return statistics.correlation(ra, rb)


markers, gene2id, factors = {}, {}, {}
sets = {}
for line in open(sys.argv[1]):
    if line.startswith(("#", "set\t")):
        continue
    mset, marker, name, length = line.rstrip("\n").split("\t")[:4]
    markers[marker] = (name, int(length))
    sets.setdefault(marker, []).append(mset)
for line in open(sys.argv[2]):
    marker, gid = line.split()
    gene2id[marker] = int(gid)
head = None
for line in open(sys.argv[3]):
    f = line.rstrip("\n").split("\t")
    if head is None:
        head = f
        continue
    row = dict(zip(head, f))
    factors[int(row["geneid"])] = row

rows = []
for marker, (name, length) in markers.items():
    row = factors[gene2id[marker]]
    rows.append((marker, name, "+".join(sets[marker]), marker_category(name), CATEGORY_RATE[marker_category(name)],
                 float(row["within_factor"]), float(row["between_factor"]), int(row["species"]), length))
print(f"{len(rows)} genes; Spearman of the category rate with the r226 within-species factor "
      f"{spearman([r[4] for r in rows], [r[5] for r in rows]):+.3f}, with the between-congener factor "
      f"{spearman([r[4] for r in rows], [r[6] for r in rows]):+.3f}")
print("\ncategory      genes  rate  r226 within factor median (quartiles)   between median")
for cat in CATEGORY_RATE:
    sel = [r for r in rows if r[3] == cat]
    w = sorted(r[5] for r in sel)
    q = statistics.quantiles(w, n=4) if len(w) > 1 else [w[0]] * 3
    print(f"{cat:12s} {len(sel):6d} {CATEGORY_RATE[cat]:5.1f}  {statistics.median(w):.2f} ({q[0]:.2f}-{q[2]:.2f})"
          f"{'':20s}{statistics.median(r[6] for r in sel):.2f}")
print("\nfactor spread: real within", f"{min(r[5] for r in rows):.2f}-{max(r[5] for r in rows):.2f}",
      "(CV", f"{statistics.pstdev(r[5] for r in rows) / statistics.mean(r[5] for r in rows):.2f})")
print("\nper set: species with the gene in r226 (gene_congeners 'species'), median")
for s in ("bac120", "ar53", "bac120+ar53"):
    sel = [r for r in rows if r[2] == s]
    print(f"  {s:12s} {len(sel):4d} genes, median species {statistics.median(r[7] for r in sel):.0f}")
print("\ngenes the categories place worst (|log2(rate / factor)| largest):")
import math  # noqa: E402
worst = sorted(rows, key=lambda r: -abs(math.log2(r[4] / r[5])))[:15]
print("marker        name                 set           category     rate  r226 factor")
for r in worst:
    print(f"{r[0]:13s} {r[1][:20]:20s} {r[2]:13s} {r[3]:12s} {r[4]:4.1f}  {r[5]:.2f}")
with open("rates_vs_r226.tsv", "w") as fh:
    fh.write("marker\tname\tset\tcategory\tcategory_rate\tr226_within_factor\tr226_between_factor\tr226_species\thmm_length\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
