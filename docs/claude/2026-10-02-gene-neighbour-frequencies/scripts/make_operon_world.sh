#!/bin/bash
# The tuning world's design (900 species, 8% archaea, 3 genomes each, strains 0.2-2%, congeners 1.5-6%, 15% of
# species unknown to every database) with operon-like marker clusters, each family breaking them up (--operons,
# --operon_breaks 0.25): gene order differs between clades. For testing adjacency evidence in the presence model.
set -e
T=$HOME/opw; R=/mnt/c/Users/hildebra/Documents/locDev/protal
mkdir -p $T
python3 $R/scripts/mini_db/gtdb_like_lineages.py --species 900 --archaea 0.08 --seed 1 > $T/lineages.txt
/usr/bin/time -f "release %e s" python3 $R/scripts/mini_db/simulate_gtdb_release.py --outdir $T/world --lineages $T/lineages.txt \
    --seed 21 --genomes_per_species 3 --strain_divergence 0.002-0.02 --species_divergence 0.015-0.06 \
    --operons --operon_breaks 0.25 2> $T/world.log
tail -2 $T/world.log
python3 $R/docs/claude/2026-09-29-model-training-tuning/make_release_p.py $T/world $T/release_p $T/unknown.txt 0.15 7
du -sh $T/world
echo WORLD_DONE
