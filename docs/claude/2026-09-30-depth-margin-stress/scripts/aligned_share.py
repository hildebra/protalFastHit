#!/usr/bin/env python3
"""aligned_share.py - aligned mates (primary records, db_full) per simulated mate, by the genome's distance
from the reference: what a strain loses before any depth rule applies. Reads ~/stress."""
import collections, csv, os, re, subprocess, statistics
B = os.path.expanduser("~/stress")
design = collections.defaultdict(dict)
for r in csv.DictReader(open(f"{B}/design/design.tsv"), delimiter="\t"):
    design[r["sample"]][r["genome"]] = r
pairs = {}
for rl in (100, 150):
    for r in csv.DictReader(open(f"{B}/design/rl{rl}.manifest.tsv"), delimiter="\t"):
        pairs[(r["sample"], r["genome"])] = int(r["read_pairs"])
ratio = collections.defaultdict(list)
for sample in sorted(design):
    counts = collections.Counter()
    p = subprocess.Popen(["zstd", "-dc", f"{B}/align/full/{sample}/{sample}.sam.zst"], stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        if line.startswith("@"): continue
        f = line.split("\t", 3)
        if int(f[1]) & 0x904: continue
        counts[re.match(r"(GC[AF]_\d{9}\.\d+)", f[0]).group(1)] += 1
    for g, r in design[sample].items():
        d = float(r["distance"])
        b = "reference" if d == 0 else "<2%" if d < 0.02 else "2-4%" if d < 0.04 else "4-5%" if d < 0.05 else ">=5%"
        ratio[(sample[:4], b)].append(counts[g] / (2 * pairs[(sample, g)]))
print("aligned mates per simulated mate (primary records, db_full), median over genomes")
for key in sorted(ratio):
    print(key, f"{statistics.median(ratio[key]):.3f}", f"n={len(ratio[key])}")
