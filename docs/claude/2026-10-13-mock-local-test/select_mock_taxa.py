"""Choose the GTDB r226 genomes of a reduced database for a local test of the mock benchmark (scripts/mock_benchmark):
every species of the truth tables (scripts/mock_benchmark/truth/*.tsv), up to CONGENERS more species of each of their
genera as competitors (drawn by a seeded hash), and up to STRAINS other genomes per species (CheckM2 >= 90% complete,
<= 5% contaminated; isolates first, then by the hash). Writes OUT/bac120_taxa.tsv and OUT/ar53_taxa.tsv in the format
of docs/claude/2026-10-10-real-ancestry/extract_subset.py (accession, role, species, genus, family).

usage: select_mock_taxa.py <GTDB release dir> <truth dir> <out dir> [congeners] [strains] [seed]
"""
import glob
import gzip
import hashlib
import os
import sys
from collections import defaultdict

REL, TRUTH, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
CONGENERS = int(sys.argv[4]) if len(sys.argv) > 4 else 3
STRAINS = int(sys.argv[5]) if len(sys.argv) > 5 else 5
SEED = sys.argv[6] if len(sys.argv) > 6 else "1"


def h(*parts):
    return hashlib.sha1(("|".join(map(str, parts)) + "|" + SEED).encode()).hexdigest()


wanted = set()
for path in glob.glob(os.path.join(TRUTH, "*.tsv")):
    with open(path) as f:
        rows = [line.rstrip("\n").split("\t") for line in f if not line.startswith("#")]
    wanted |= {r[0] for r in rows[1:] if r[0].startswith("s__")}

genomes = defaultdict(list)  # species -> [(accession, is_rep, isolate)]
lineage, marker_set = {}, {}
for SET in ("bac120", "ar53"):
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
                isolate = r[col["ncbi_genome_category"]] in ("", "none", "derived from isolate")
                genomes[sp].append((r[col["accession"]], rep, isolate))
                lineage[sp] = tax
                marker_set[sp] = SET

missing = sorted(wanted - set(lineage))
if missing:
    sys.exit(f"truth species not in r226: {missing}")
by_genus = defaultdict(list)
for sp, tax in lineage.items():
    if any(rep for _, rep, _ in genomes[sp]):
        by_genus[tax[5]].append(sp)
chosen = set(wanted)
for genus in sorted({lineage[sp][5] for sp in wanted}):
    others = sorted((s for s in by_genus[genus] if s not in wanted), key=lambda s: h("congener", s))
    chosen |= set(others[:CONGENERS])

rows = {"bac120": [], "ar53": []}
for sp in sorted(chosen):
    tax = lineage[sp]
    isolate = {a: i for a, _, i in genomes[sp]}
    strains = sorted((a for a, rep, _ in genomes[sp] if not rep), key=lambda a: (not isolate[a], h("strain", a)))
    for a in (a for a, rep, _ in genomes[sp] if rep):
        rows[marker_set[sp]].append((a, "representative", sp, tax[5], tax[4]))
    for a in strains[:STRAINS]:
        rows[marker_set[sp]].append((a, "strain", sp, tax[5], tax[4]))

os.makedirs(OUT, exist_ok=True)
for SET, rs in rows.items():
    with open(os.path.join(OUT, f"{SET}_taxa.tsv"), "w") as out:
        out.write("accession\trole\tspecies\tgenus\tfamily\n")
        for r in rs:
            out.write("\t".join(r) + "\n")
    print(f"{SET}: {len({r[2] for r in rs})} species ({len({r[2] for r in rs} & wanted)} of the truth), "
          f"{len({r[3] for r in rs})} genera, {sum(r[1] == 'strain' for r in rs)} strains, {len(rs)} genomes")
print(f"truth species {len(wanted)}, chosen {len(chosen)} (congeners up to {CONGENERS} per genus, strains up to "
      f"{STRAINS} per species, seed {SEED})")
