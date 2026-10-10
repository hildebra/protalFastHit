"""Choose real GTDB r226 taxa for testing the ancestry sites, strain alleles and column weights on real marker genes:
FAMILIES bacterial families, GENERA genera each, every species of those genera (so that a species' congeners are all
its real congeners), and up to MAX_STRAINS other genomes per species (CheckM2 >= 90% complete, <= 5% contaminated;
chosen by a hash of the accession, isolates first).

A genus qualifies with 4 to 12 species of which at least 2 have 3 or more genomes (the representative and two
strains); a family with GENERA or more qualifying genera. Families and genera are drawn by a seeded hash.

usage: select_taxa.py <GTDB release dir (bac120_metadata_r226.tsv.gz)> <out taxa.tsv> [seed] [marker set: bac120, ar53]
"""
import gzip
import hashlib
import sys
from collections import defaultdict

REL, OUT = sys.argv[1], sys.argv[2]
SEED = sys.argv[3] if len(sys.argv) > 3 else "1"
SET = sys.argv[4] if len(sys.argv) > 4 else "bac120"
FAMILIES, GENERA = 5, 4
MIN_SPECIES, MAX_SPECIES, MIN_STRAIN_SPECIES, MIN_GENOMES = 4, 12, 2, 3
MAX_STRAINS = 20


def h(*parts):
    return hashlib.sha1(("|".join(map(str, parts)) + "|" + SEED).encode()).hexdigest()


genomes = defaultdict(list)  # species -> [(accession, is_rep, isolate)]
lineage = {}
with gzip.open(f"{REL}/{SET}_metadata_r226.tsv.gz", "rt") as f:
    head = f.readline().rstrip("\n").split("\t")
    col = {c: i for i, c in enumerate(head)}
    for line in f:
        r = line.rstrip("\n").split("\t")
        tax = r[col["gtdb_taxonomy"]].split(";")
        sp = tax[6]
        rep = r[col["gtdb_representative"]] == "t"
        try:
            ok = float(r[col["checkm2_completeness"]]) >= 90 and float(r[col["checkm2_contamination"]]) <= 5
        except ValueError:
            ok = False
        if rep or ok:
            genomes[sp].append((r[col["accession"]], rep, r[col["ncbi_genome_category"]] in ("", "none", "derived from isolate")))
            lineage[sp] = tax

by_genus = defaultdict(list)
for sp, tax in lineage.items():
    if any(rep for _, rep, _ in genomes[sp]):
        by_genus[(tax[4], tax[5])].append(sp)


def qualifies(species):
    rich = sum(len(genomes[s]) >= MIN_GENOMES for s in species)
    return MIN_SPECIES <= len(species) <= MAX_SPECIES and rich >= MIN_STRAIN_SPECIES


families = defaultdict(list)
for (fam, gen), species in by_genus.items():
    if qualifies(species) and not gen.endswith("__"):
        families[fam].append(gen)
eligible = sorted((f for f, g in families.items() if len(g) >= GENERA), key=lambda f: h("family", f))
chosen = eligible[:FAMILIES]

rows = []
for fam in chosen:
    for gen in sorted(families[fam], key=lambda g: h("genus", g))[:GENERA]:
        for sp in sorted(by_genus[(fam, gen)]):
            reps = [a for a, rep, _ in genomes[sp] if rep]
            strains = sorted((a for a, rep, iso in genomes[sp] if not rep),
                             key=lambda a: (not dict((x, i) for x, _, i in genomes[sp])[a], h("strain", a)))
            for a in reps:
                rows.append((a, "representative", sp, gen, fam))
            for a in strains[:MAX_STRAINS]:
                rows.append((a, "strain", sp, gen, fam))

with open(OUT, "w") as out:
    out.write("accession\trole\tspecies\tgenus\tfamily\n")
    for r in rows:
        out.write("\t".join(r) + "\n")
print(f"{len(eligible)} families qualify; chosen: {', '.join(chosen)}")
for fam in chosen:
    fr = [r for r in rows if r[4] == fam]
    gens = sorted({r[3] for r in fr})
    print(f"  {fam}: {len(gens)} genera, {len({r[2] for r in fr})} species, "
          f"{sum(r[1] == 'strain' for r in fr)} strains ({', '.join(g[3:] for g in gens)})")
print(f"total: {len({r[3] for r in rows})} genera, {len({r[2] for r in rows})} species, "
      f"{sum(r[1] == 'strain' for r in rows)} strains, {len(rows)} genomes")
