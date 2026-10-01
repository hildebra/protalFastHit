#!/usr/bin/env python3
"""congener_design.py WORLD_DIR - held-out species and samples of the congener world, for simulate_metagenomes
--from_manifest.

The world (congener_world.sh): 160 species, 120 of them in 12 genera of 10; 5 genomes per species, the
representative (the database's reference) and 4 other strains, each genome 0.2-2% from its species' ancestor
(a strain 0.4-4% from the reference); each species 1-5% from its genus' ancestor at a gene of factor 1
(congeners 2-10% apart), genes of different conservation (--gene_rates categories).

Held out (heldout.txt, left out of db_missing): 4 of the 10 species of every large genus.

Samples: for each read setup (2x100 HS20 300+-40, 2x150 HS25 350+-50), 6 samples, each with every species of
4 large genera (40 species, held-out ones included) and 20 species of the small genera. Each species is, in
turn, the representative alone, one other strain alone (the strains in turn), the representative with the
species' most distant strain as a 10-30% minor (close_major), or the reverse (far_major). A species' depth is
lognormal (median 10x, sigma 1.0, 1.5-120x); in samples 4-6 of each setup the held-out species have 5 times
that, so that their reads outnumber those of their congeners in the database.
Writes design/rl<length>.manifest.tsv and design/design.tsv: sample, species, kind, genome, role (alone, major,
minor), depth, distance of the genome to the reference, held out, minor fraction.
"""
import collections
import csv
import math
import os
import random
import sys

B = sys.argv[1]
SIM = os.path.join(B, "gtdb", "simulation")
SETUPS = [(100, "s100"), (150, "s150")]
SAMPLES, LARGE_GENERA, SMALL_SPECIES = 6, 4, 20
KINDS = ["alone_ref", "alone_strain", "close_major", "far_major"]
MINOR = [0.1, 0.2, 0.3]
rng = random.Random(2027)

genomes = {}
with open(os.path.join(SIM, "genomes.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        row["species"] = row["gtdb_taxonomy"].split(";s__")[-1]
        genomes[row["accession"]] = row
div = {}
with open(os.path.join(SIM, "divergence.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        div[row["accession"]] = float(row["strain_divergence"])
by_species = collections.defaultdict(list)
for acc in sorted(genomes):
    by_species[genomes[acc]["species"]].append(acc)
rep = {s: next(a for a in accs if a.startswith("GCF_")) for s, accs in by_species.items()}
distance = {a: 0.0 if a == rep[genomes[a]["species"]] else div[a] + div[rep[genomes[a]["species"]]] for a in genomes}
strains = {s: sorted((a for a in accs if a != rep[s]), key=lambda a: distance[a]) for s, accs in by_species.items()}
species = sorted(by_species)
genus = {s: genomes[rep[s]]["gtdb_taxonomy"].split(";")[5] for s in species}
by_genus = collections.defaultdict(list)
for s in species:
    by_genus[genus[s]].append(s)
large = sorted(g for g, members in by_genus.items() if len(members) >= 10)
small_species = sorted(s for s in species if genus[s] not in large)

heldout = set()
for g in large:
    heldout.update(rng.sample(sorted(by_genus[g]), 4))
with open(os.path.join(B, "heldout.txt"), "w") as fh:
    fh.writelines(f"s__{s}\n" for s in sorted(heldout))

COLS = ["sample", "genome", "species", "taxonomy", "genome_length", "read_pairs", "vertical_coverage",
        "relative_abundance", "fasta_path", "art_seed"]
os.makedirs(os.path.join(B, "design"), exist_ok=True)
design = []
seed = 7000
for length, prefix in SETUPS:
    rows = []
    for n in range(1, SAMPLES + 1):
        sample = f"{prefix}_{n}"
        chosen = [s for g in rng.sample(large, LARGE_GENERA) for s in sorted(by_genus[g])]
        chosen += rng.sample(small_species, SMALL_SPECIES)
        rng.shuffle(chosen)
        dominant_missing = n > SAMPLES // 2
        for i, s in enumerate(chosen):
            kind = KINDS[i % len(KINDS)]
            depth = min(120.0, max(1.5, math.exp(rng.gauss(math.log(10), 1.0))))
            if dominant_missing and s in heldout:
                depth *= 5
            minor = MINOR[(i // len(KINDS) + n) % len(MINOR)]
            far = strains[s][-1]
            if kind == "alone_ref":
                parts = [(rep[s], "alone", depth)]
            elif kind == "alone_strain":
                parts = [(strains[s][(i // len(KINDS) + n) % len(strains[s])], "alone", depth)]
            elif kind == "close_major":
                parts = [(rep[s], "major", depth * (1 - minor)), (far, "minor", depth * minor)]
            else:
                parts = [(far, "major", depth * (1 - minor)), (rep[s], "minor", depth * minor)]
            for acc, role, d in parts:
                g = genomes[acc]
                size = int(g["genome_length"])
                pairs = max(1, round(d * size / (2 * length)))
                seed += 1
                rows.append(dict(sample=sample, genome=acc, species=s, taxonomy=g["gtdb_taxonomy"], genome_length=size,
                                 read_pairs=pairs, vertical_coverage=f"{pairs * 2 * length / size:.4f}",
                                 relative_abundance="", fasta_path=g["fasta_path"], art_seed=seed))
                design.append(dict(sample=sample, species=s, kind=kind, genome=acc, role=role,
                                   depth=f"{pairs * 2 * length / size:.4f}", distance=f"{distance[acc]:.4f}",
                                   heldout=int(s in heldout), minor_fraction=minor if role != "alone" else ""))
    total = collections.Counter()
    for r in rows:
        total[r["sample"]] += float(r["vertical_coverage"])
    for r in rows:
        r["relative_abundance"] = f"{float(r['vertical_coverage']) / total[r['sample']]:.6g}"
    with open(os.path.join(B, "design", f"rl{length}.manifest.tsv"), "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        fh.writelines("\t".join(str(r[c]) for c in COLS) + "\n" for r in rows)
with open(os.path.join(B, "design", "design.tsv"), "w") as fh:
    cols = list(design[0])
    fh.write("\t".join(cols) + "\n")
    fh.writelines("\t".join(str(d[c]) for c in cols) + "\n" for d in design)
print(f"{len(species)} species in {len(by_genus)} genera ({len(large)} of 10), {len(heldout)} held out", file=sys.stderr)
