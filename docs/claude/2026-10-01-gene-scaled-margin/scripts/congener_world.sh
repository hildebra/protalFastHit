#!/usr/bin/env bash
# The congener world in $WORLD (default ~/stress3): 160 species, 120 of them in 12 genera of 10
# (congener_lineages.py), 5 genomes each (simulate_gtdb_release.py: every genome 0.2-2% from its species'
# ancestor, a strain 0.4-4% from the representative; each species 1-5% from its genus' ancestor, congeners
# 2-10% apart at a gene of factor 1; genes of different conservation, --gene_rates categories); the design
# (congener_design.py). The databases and samples are made by run_rules.sh.
set -euo pipefail
B=${WORLD:-$HOME/stress3}
SRC=${SRC:-$HOME/fix-build/src}
HERE=$(cd $(dirname $0) && pwd)
S=$SRC/scripts/mini_db
mkdir -p $B/logs
[ -s $B/lineages.txt ] || python3 $HERE/congener_lineages.py 3 > $B/lineages.txt
if [ ! -f $B/gtdb/simulation/genomes.tsv ]; then
  python3 $S/simulate_gtdb_release.py --outdir $B/gtdb --lineages $B/lineages.txt --genomes_per_species 5 \
    --genome_length 200000 --strain_divergence 0.002-0.02 --species_divergence 0.01-0.05 --gene_rates categories \
    --seed 13 > $B/logs/release.log 2>&1
fi
python3 $HERE/congener_design.py $B
echo "world ready in $B"
