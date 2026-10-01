#!/usr/bin/env bash
# Each 0.7 version's own pipeline on the benchmark world's release (765 species; the 135 unknown ones are in
# no database): build_gtdb_database.py of 0.7.0 (4b21427) with protal 0.7.0, and of 0.7.1 (1c11a00) with
# protal 0.7.1, one after the other, the same design as the V2 run of 2026-09-30 (reduced sizes), seed 1.
# Each writes its finished database (protal_db: all species, its four trained models; 0.7.1's also
# gene_neighbours.tsv and gene_conservation.tsv) and its training database (training_db: held-out species
# and clades left out). simulate_metagenomes of 0.7.1 for both (the samples the models are trained on).
set -uo pipefail
B=${BENCH:-$HOME/bench071}
W=$B/world
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
PBSIM=${PBSIM:-$HOME/micromamba/envs/protal-db-build/bin/pbsim}
T=${T:-6}
for v in ${VERSIONS:-0.7.0 0.7.1}; do
  out=$B/V${v//./}
  [ -f $out/protal_db/database.protal ] && [ -s $out/model_logs/summary.txt ] && { echo "$out done"; continue; }
  echo "$(date +%T) pipeline of $v into $out"
  /usr/bin/time -v -o $B/logs/pipeline_$v.time $PY $B/src/$v/scripts/build_gtdb_database.py --gtdb $W/release_p \
    --outdir $out --protal $B/bin/protal-$v --simulator $B/bin/simulate_metagenomes -t $T --seed 1 \
    --extra-genomes $W/full/simulation/genomes_nonreps --samples 4 --read-pairs 1000,20000,200000 \
    --test-samples 2 --test-read-pairs 500,10000,500000 --long-read-bases 300000,6000000,60000000 \
    --test-long-read-bases 150000,3000000,90000000 --pbsim $PBSIM > $B/logs/pipeline_$v.log 2>&1 ||
    { echo "pipeline of $v failed, see $B/logs/pipeline_$v.log"; tail -20 $B/logs/pipeline_$v.log; exit 1; }
  grep -E "Elapsed|Maximum resident" $B/logs/pipeline_$v.time
  cat $out/model_logs/summary.txt
done
