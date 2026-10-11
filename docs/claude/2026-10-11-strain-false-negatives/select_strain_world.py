"""Choose a strain-rich world from GTDB r226 for the strain false-negative study: FOCAL bacterial species per stratum of
GTDB's mean intra-species ANI (tight >= 99.0, mid 97.5-99.0, wide < 97.5), each from its own genus of 2 to 30 species,
with at least MIN_OTHERS genomes besides the representative (CheckM2 >= 90% complete, <= 5% contaminated). Of each
focal species up to 2 * HALF other genomes, ordered by a seeded hash and dealt alternately to the allele half (may go
into a database's full reference) and the test half (never in a database; reads are simulated from them). Every
other species of the focal genera comes with its representative (the congeners).

Writes taxa.tsv (accession, role, species, genus, family, stratum, rank) for extract_subset.py, whose roles are
representative and strain; the role column says allele/test in a second file, halves.tsv.

usage: select_strain_world.py <GTDB release dir> <out dir> [seed]
"""
import gzip
import hashlib
import os
import sys
from collections import defaultdict

REL, OUT = sys.argv[1], sys.argv[2]
SEED = sys.argv[3] if len(sys.argv) > 3 else "1"
FOCAL, HALF, MIN_OTHERS = 10, 20, 24
GENUS_MIN, GENUS_MAX = 2, 30


def h(*parts):
    return hashlib.sha1(("|".join(map(str, parts)) + "|" + SEED).encode()).hexdigest()


genomes = defaultdict(list)  # species -> [(accession, is_rep)]
lineage = {}
with gzip.open(f"{REL}/bac120_metadata_r226.tsv.gz", "rt") as f:
    col = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for line in f:
        r = line.rstrip("\n").split("\t")
        tax = r[col["gtdb_taxonomy"]].split(";")
        rep = r[col["gtdb_representative"]] == "t"
        try:
            ok = float(r[col["checkm2_completeness"]]) >= 90 and float(r[col["checkm2_contamination"]]) <= 5
        except ValueError:
            ok = False
        if rep or ok:
            genomes[tax[6]].append((r[col["accession"]], rep))
            lineage[tax[6]] = tax

mean_ani = {}
with open(f"{REL}/auxillary_files/sp_clusters_r226.tsv") as f:
    col = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for line in f:
        r = line.rstrip("\n").split("\t")
        try:
            mean_ani[r[col["GTDB species"]]] = float(r[col["Mean intra-species ANI"]])
        except ValueError:
            pass

by_genus = defaultdict(list)
for sp, tax in lineage.items():
    if any(rep for _, rep in genomes[sp]) and not tax[5].endswith("__"):
        by_genus[tax[5]].append(sp)


def stratum(sp):
    a = mean_ani.get(sp)
    if a is None:
        return None
    return "tight" if a >= 99.0 else "mid" if a >= 97.5 else "wide"


candidates = defaultdict(list)
for gen, species in by_genus.items():
    if not GENUS_MIN <= len(species) <= GENUS_MAX:
        continue
    for sp in species:
        others = [a for a, rep in genomes[sp] if not rep]
        if len(others) >= MIN_OTHERS and stratum(sp):
            candidates[stratum(sp)].append(sp)

focal, used_genera = [], set()
for st in ("tight", "mid", "wide"):
    n = 0
    for sp in sorted(candidates[st], key=lambda s: h("species", s)):
        if lineage[sp][5] in used_genera:
            continue
        focal.append((sp, st))
        used_genera.add(lineage[sp][5])
        n += 1
        if n == FOCAL:
            break

os.makedirs(OUT, exist_ok=True)
rows, halves = [], []
for gen in sorted(used_genera):
    for sp in sorted(by_genus[gen]):
        st = dict(focal).get(sp, "congener")
        tax = lineage[sp]
        for a, rep in genomes[sp]:
            if rep:
                rows.append((a, "representative", sp, gen, tax[4], st, "0"))
        if st == "congener":
            continue
        others = sorted((a for a, rep in genomes[sp] if not rep), key=lambda a: h("strain", a))[:2 * HALF]
        for i, a in enumerate(others):
            rows.append((a, "strain", sp, gen, tax[4], st, str(i // 2 + 1)))
            halves.append((a, "allele" if i % 2 == 0 else "test", sp, st, str(i // 2 + 1)))

with open(f"{OUT}/taxa.tsv", "w") as o:
    o.write("accession\trole\tspecies\tgenus\tfamily\tstratum\trank\n")
    o.writelines("\t".join(r) + "\n" for r in rows)
with open(f"{OUT}/halves.tsv", "w") as o:
    o.write("accession\thalf\tspecies\tstratum\trank\n")
    o.writelines("\t".join(r) + "\n" for r in halves)
for sp, st in focal:
    n = sum(1 for r in halves if r[2] == sp)
    print(f"{st}\t{sp}\tmean ANI {mean_ani[sp]:.2f}\t{len(by_genus[lineage[sp][5]])} congeners+1\t{n} strains\t"
          f"{lineage[sp][1]}")
print(f"total: {len(used_genera)} genera, {len({r[2] for r in rows})} species, {len(halves)} strains, {len(rows)} genomes")
