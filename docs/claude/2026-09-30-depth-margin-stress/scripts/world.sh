#!/usr/bin/env bash
# The stress world: 120 GTDB-like species with 5 genomes each, every genome 0.2-3% from its species'
# ancestor (a strain 0.4-6% from the representative), congeneric species 3-12% apart at the markers; the
# design (design.py); and the database built by 0.7 twice, with every species (db_full) and without the
# held-out quarter (db_missing), at zstd level 3. $RELEASE_ARGS go to simulate_gtdb_release.py (the second world:
# --gene_rates categories).
set -euo pipefail
B=${STRESS:-$HOME/stress}
SRC7=${SRC7:-$HOME/fix-build/src}
P7=${P7:-$HOME/fix-build/bin/protal}
T=${T:-6}
HERE=$(cd $(dirname $0) && pwd)
S=$SRC7/scripts/mini_db
mkdir -p $B/logs
[ -s $B/lineages.txt ] || python3 $S/gtdb_like_lineages.py --species 120 --archaea 0.1 --seed 5 > $B/lineages.txt
if [ ! -f $B/gtdb/simulation/genomes.tsv ]; then
  python3 $S/simulate_gtdb_release.py --outdir $B/gtdb --lineages $B/lineages.txt --genomes_per_species 5 \
    --genome_length 200000 --strain_divergence 0.002-0.03 --species_divergence 0.015-0.06 --seed 11 ${RELEASE_ARGS:-} > $B/logs/release.log 2>&1
fi
python3 $HERE/design.py $B
for db in full missing; do
  [ -f $B/db_$db/database.protal ] && continue
  extra=$([ $db = missing ] && echo "--exclude_species $B/heldout.txt" || true)
  python3 $S/gtdb_to_protal_db.py --gtdb $B/gtdb --outdir $B/db_$db -t $T $extra > $B/logs/convert_$db.log 2>&1
  $P7 --build --no_profile -t $T --db $B/db_$db --reference $B/db_$db/reference.fna \
    --full_reference $B/db_$db/full_reference.fna --compress_level 3 > $B/logs/build_$db.log 2>&1
done
echo "world ready in $B"
