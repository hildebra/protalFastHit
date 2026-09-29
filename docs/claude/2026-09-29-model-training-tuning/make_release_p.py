#!/usr/bin/env python3
"""make_release_p.py WORLD OUT UNKNOWN.txt FRACTION SEED: a release like WORLD without a random FRACTION of its
species (per domain; written to UNKNOWN.txt): those species exist in no database, like the species GTDB lacks.
Taxonomy and metadata are filtered; the rest is linked, since the converter and the genome table keep only
genomes in the taxonomy."""
import collections, gzip, os, random, sys
world, out, unknown_path, fraction, seed = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4]), int(sys.argv[5])
rel = "r226"
lineage = {}
for mset in ("bac120", "ar53"):
    for line in open(os.path.join(world, f"{mset}_taxonomy_{rel}.tsv")):
        acc, lin = line.rstrip("\n").split("\t")[:2]
        lineage[acc] = lin
by_domain = collections.defaultdict(set)
for lin in lineage.values():
    by_domain[lin.split(";")[0]].add(lin.split(";")[-1])
rng = random.Random(seed)
unknown = set()
for d in sorted(by_domain):
    unknown |= set(rng.sample(sorted(by_domain[d]), round(fraction * len(by_domain[d]))))
os.makedirs(out, exist_ok=True)
for name in os.listdir(world):
    src, dst = os.path.join(world, name), os.path.join(out, name)
    if os.path.lexists(dst):
        os.remove(dst)
    if "_taxonomy_" in name:
        with open(dst, "w") as fh:
            fh.writelines(l for l in open(src) if l.split("\t")[1].strip().split(";")[-1] not in unknown)
    elif "_metadata_" in name:
        with gzip.open(src, "rt") as fin, gzip.open(dst, "wt") as fout:
            header = fin.readline()
            fout.write(header)
            col = header.rstrip("\n").split("\t").index("gtdb_taxonomy")
            fout.writelines(l for l in fin if l.rstrip("\n").split("\t")[col].split(";")[-1] not in unknown)
    elif name != "simulation":
        os.symlink(os.path.abspath(src), dst)
open(unknown_path, "w").writelines(s + "\n" for s in sorted(unknown))
print(f"{len(unknown)} of {sum(len(v) for v in by_domain.values())} species unknown to every database -> {unknown_path}")
