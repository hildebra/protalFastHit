#!/bin/bash
# The GTDB-like tuning world: 900 species (8% archaea), 3 genomes each, strains 0.4-4% from their
# representative, congeneric species 3-12% apart at the markers; 15% of species unknown to every database.
set -e
T=$HOME/tune; R=/mnt/c/Users/hildebra/Documents/locDev/protal
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/fbe611b7-59a8-4a3b-b159-58842b2426b7/scratchpad/tune
mkdir -p $T
python3 $R/scripts/mini_db/gtdb_like_lineages.py --species 900 --archaea 0.08 --seed 1 > $T/lineages.txt
/usr/bin/time -f "release %e s" python3 $R/scripts/mini_db/simulate_gtdb_release.py --outdir $T/world --lineages $T/lineages.txt \
    --seed 21 --genomes_per_species 3 --strain_divergence 0.002-0.02 --species_divergence 0.015-0.06 2> $T/world.log
tail -2 $T/world.log
python3 $S/make_release_p.py $T/world $T/release_p $T/unknown.txt 0.15 7
du -sh $T/world
