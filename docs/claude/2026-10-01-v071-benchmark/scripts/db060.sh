#!/usr/bin/env bash
# 0.6.0a's databases, built by 0.6.0a from the release as 0.7.1's converter writes it (the converter of the
# pipelines; full_reference.fna decompressed, which 0.6.0a needs plain): db060_full (765 species) and
# db060_missing (without the species and clades 0.7.1's pipeline held out of its training database,
# V071/heldout_species.txt), each with the model 0.6.0a ships (scripts/random_forest.xml as model.xml).
set -uo pipefail
B=${BENCH:-$HOME/bench071}
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
CONV=$B/src/0.7.1/scripts/mini_db/gtdb_to_protal_db.py
for kind in full missing; do
  db=$B/db060_$kind
  [ -s $db/index.prx ] && continue
  rm -rf $db
  extra=()
  [ $kind = missing ] && extra=(--exclude_species $B/V071/heldout_species.txt)
  $PY $CONV --gtdb $B/world/release_p --outdir $db -t 6 "${extra[@]}" > $B/logs/convert060_$kind.log 2>&1 || { echo "convert $kind failed"; exit 1; }
  [ -f $db/full_reference.fna.zst ] && zstd -dq --rm $db/full_reference.fna.zst
  cp $B/src/0.6.0a/scripts/random_forest.xml $db/model.xml
  /usr/bin/time -v -o $B/logs/build060_$kind.time $B/bin/protal-0.6.0a --build --no_profile -t 6 --db $db \
    --reference $db/reference.fna --full_reference $db/full_reference.fna > $B/logs/build060_$kind.log 2>&1 || { echo "build $kind failed"; exit 1; }
  grep -E "Elapsed|Maximum resident" $B/logs/build060_$kind.time
done
