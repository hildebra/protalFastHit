#!/bin/bash
# A dense world: GENERA genera of SPECIES_PER_GENUS species each (default 4 x 60), so that a marker gene has up to
# 60 relatives 1.5-6% apart in its genus, as in the big GTDB genera, where the 900-species world has mostly 1-3
# species per genus. Writes the release into $PERF_DIR/dense/world, then the protal database and reads:
#   dense_world.sh [GENERA] [SPECIES_PER_GENUS]
set -e
source "$(dirname "$0")/env.sh"
genera=${1:-4}; per=${2:-60}
D=$PERF_DIR/dense; rm -rf $D; mkdir -p $D
awk -v g=$genera -v n=$per 'BEGIN { for (i = 1; i <= g; i++) for (j = 1; j <= n; j++) printf "d__Bacteria;p__Phylum%d;c__Class%d;o__Order%d;f__Family%d;g__Genus%d;s__Genus%d sp%03d\n", i, i, i, i, i, i, j }' > $D/lineages.txt
python3 $PROTAL_SRC/scripts/mini_db/simulate_gtdb_release.py --outdir $D/world --lineages $D/lineages.txt \
  --strain_divergence 0.002-0.02 --species_divergence 0.015-0.06 > $D/world.log 2>&1
tail -2 $D/world.log
