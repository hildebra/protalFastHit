#!/usr/bin/env bash
# The benchmark world: 900 GTDB-like species (12% archaea), 3 genomes each, strains 0.4-4% from their
# representative, congeneric species 3-12% apart at a gene of factor 1, genes of different conservation
# (--gene_rates categories) and markers in operon-like clusters (--operons, new in 0.7.1's simulator); 15% of
# the species are unknown to every database (release_p, by the tuning study's make_release_p.py). Scripts of
# 0.7.1 (the archive in $B/src/0.7.1).
set -euo pipefail
B=${BENCH:-$HOME/bench071}
S=$B/src/0.7.1/scripts/mini_db
TUNE=$B/src/0.7.1/docs/claude/2026-09-29-model-training-tuning
W=$B/world
mkdir -p $W $B/logs
[ -s $W/lineages.txt ] || python3 $S/gtdb_like_lineages.py --species 900 --archaea 0.12 --seed 3 > $W/lineages.txt
if [ ! -f $W/full/simulation/genomes.tsv ]; then
  /usr/bin/time -f "release %e s" python3 $S/simulate_gtdb_release.py --outdir $W/full --lineages $W/lineages.txt \
    --seed 31 --genomes_per_species 3 --strain_divergence 0.002-0.02 --species_divergence 0.015-0.06 \
    --gene_rates categories --operons > $B/logs/release.log 2>&1
fi
[ -d $W/release_p ] || python3 $TUNE/make_release_p.py $W/full $W/release_p $W/unknown.txt 0.15 7
tail -2 $B/logs/release.log
du -sh $W/full
