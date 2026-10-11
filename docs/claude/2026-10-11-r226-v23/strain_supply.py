"""What download_gtdb.py would pick from a GTDB release's metadata (its own pick(), no download) and the strain mix a
build would simulate from it: species with real strains, how often a simulated species is a real strain, an in-silico
strain or its representative, with the build's allele split (--allele-genome-share: each non-representative genome by
hash to the alleles, so k of a species' P strains stay simulated, k ~ Binomial(P, share)) and in-silico strains of
--insilico-multi of the species with real strains and --insilico-strains of the one-genome species, the latter fixed
or chosen (auto) so that real strains are at least --min-real-strains of the simulated strains.

usage: strain_supply.py <scripts dir> <GTDB folder with *_metadata_r226.tsv.gz> [release]
"""
import argparse
import math
import sys
from types import SimpleNamespace

sys.path.insert(0, sys.argv[1])
import download_gtdb  # noqa: E402

GTDB, RELEASE = sys.argv[2], (sys.argv[3] if len(sys.argv) > 3 else "226")


def binom(n, k, p):
    return math.comb(n, k) * p ** k * (1 - p) ** (n - k)


def mix(strains_per_species, n_rep_only, allele_share=0.5, multi=0.25, insilico="auto", min_real=0.5):
    """Expected draws (per simulated species) of real strains, in-silico strains and representatives."""
    real = made = 0.0
    one_genome = float(n_rep_only)
    for P in strains_per_species:
        for k in range(P + 1):
            w = binom(P, k, 1 - allele_share)  # k strains stay simulated
            if k == 0:
                one_genome += w
                continue
            real += w * ((1 - multi) * k / (k + 1) + multi * k / (k + 2))
            made += w * multi / (k + 2)
    if insilico == "auto":
        budget = real * (1 - min_real) / min_real - made
        share = max(0.0, min(1.0, budget / (0.5 * one_genome))) if one_genome else 0.0
    else:
        share = float(insilico)
    made += 0.5 * share * one_genome
    species = len(strains_per_species) + n_rep_only
    return {"species": species, "real": real / species, "insilico": made / species,
            "rep": 1 - (real + made) / species, "real_of_strains": real / max(real + made, 1e-9), "share": share}


genomes = download_gtdb.read_metadata(GTDB, RELEASE)
rows = []
print(f"{len(genomes)} genomes in the metadata", flush=True)
for species, per_species, rep_only in ((16000, 2, 9000), (16000, 2, 19000), (30000, 2, 9000), (30000, 4, 9000),
                                       (40000, 4, 5000), (40000, 4, 9000), (40000, 4, 19000), (100000, 6, 5000), (100000, 4, 0)):
    opts = SimpleNamespace(species=species, per_species=per_species, rep_only_species=rep_only, seed=1,
                           min_completeness=90.0, max_contamination=5.0, tech_candidates=30)
    picked, chosen, rep_only_list, domain, representative, candidates = download_gtdb.pick(genomes, opts)
    counts = {}
    for g, _ in picked:
        counts[g.species] = counts.get(g.species, 0) + 1
    strain_sizes = list(counts.values())
    print(f"--species {species} --per_species {per_species} --rep_only_species {rep_only}: {len(chosen)} species with "
          f"strains (of {sum(1 for s in candidates if candidates[s])} that have candidates), {len(picked)} strains, "
          f"{len(rep_only_list)} representative-only species")
    for multi, insilico in ((0.5, 1.0), (0.25, "auto"), (0.0, "auto")):
        m = mix(strain_sizes, len(rep_only_list), multi=multi, insilico=insilico)
        print(f"    --insilico-multi {multi:g} --insilico-strains {insilico}: {m['species']} species simulated; a real "
              f"strain {100 * m['real']:.1f}%, in-silico {100 * m['insilico']:.1f}% (one-genome share {m['share']:.2f}), the "
              f"representative {100 * m['rep']:.1f}%; real strains {100 * m['real_of_strains']:.1f}% of the strains")
