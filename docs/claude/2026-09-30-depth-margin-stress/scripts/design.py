#!/usr/bin/env python3
"""design.py STRESS_DIR - the stress test's held-out species and samples, for simulate_metagenomes --from_manifest.

The world (world.sh) has 5 genomes per species: the representative (the database's reference) and 4 other
strains, each genome 0.2-3% from the species' ancestor, so a strain is 0.4-6% from the reference
(simulation/divergence.tsv: the sum of the strain's and the representative's divergence).

Held out (heldout.txt, left out of db_missing): a quarter of the species, drawn among those with a congener
that stays in the database, so that their reads land on relatives.

Samples: for each read setup (2x100 HS20 300+-40, 2x150 HS25 350+-50), 6 samples of 60 species drawn at
random, 15 of each kind, in turn:
  alone_ref     the representative alone
  alone_strain  one other strain alone (the strains in turn, so every distance occurs)
  close_major   the representative plus the species' most distant strain as a minor (10, 20 or 30%)
  far_major     the most distant strain plus the representative as a minor (10, 20 or 30%)
A species' depth is lognormal (median 12x, sigma 0.8, 2-80x), split between its genomes as the kind says.
Writes design/rl<length>.manifest.tsv and design/design.tsv: sample, species, kind, genome, role (alone,
major, minor), depth, distance of the genome to the reference, held out.
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
SAMPLES, SPECIES_PER_SAMPLE = 6, 60
KINDS = ["alone_ref", "alone_strain", "close_major", "far_major"]
MINOR = [0.1, 0.2, 0.3]
rng = random.Random(2026)

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

# Held out: a quarter of the species, each with a congener that stays.
genus = {s: genomes[rep[s]]["gtdb_taxonomy"].split(";")[5] for s in species}
by_genus = collections.defaultdict(list)
for s in species:
    by_genus[genus[s]].append(s)
heldout = set()
candidates = [s for s in species if len(by_genus[genus[s]]) >= 2]
rng.shuffle(candidates)
for s in candidates:
    if len(heldout) >= len(species) // 4:
        break
    if sum(1 for x in by_genus[genus[s]] if x not in heldout and x != s) >= 1:
        heldout.add(s)
with open(os.path.join(B, "heldout.txt"), "w") as fh:
    fh.writelines(f"s__{s}\n" for s in sorted(heldout))

COLS = ["sample", "genome", "species", "taxonomy", "genome_length", "read_pairs", "vertical_coverage",
        "relative_abundance", "fasta_path", "art_seed"]
os.makedirs(os.path.join(B, "design"), exist_ok=True)
design = []
seed = 5000
for length, prefix in SETUPS:
    rows = []
    for n in range(1, SAMPLES + 1):
        sample = f"{prefix}_{n}"
        chosen = rng.sample(species, SPECIES_PER_SAMPLE)
        for i, s in enumerate(chosen):
            kind = KINDS[i % len(KINDS)]
            depth = min(80.0, max(2.0, math.exp(rng.gauss(math.log(12), 0.8))))
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
    # relative abundance: the genome's share of the sample's cells (its depth over the sample's summed depth)
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
far = [distance[strains[s][-1]] for s in species]
print(f"{len(species)} species, {len(heldout)} held out; most distant strain per species "
      f"{min(far):.3f}-{max(far):.3f} (median {sorted(far)[len(far) // 2]:.3f}) from the reference", file=sys.stderr)
