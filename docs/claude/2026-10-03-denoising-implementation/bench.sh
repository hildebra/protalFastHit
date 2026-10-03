#!/usr/bin/env bash
# Two build-and-train pipelines on the benchmark world (765 species, 135 unknown; ~/bench071/world) with the denoise
# branch's protal and scripts: "new" with the new defaults (congener groups 0.25:2-5, the relatives features, calls at a
# target share of false calls), "base" with the design, features and calls before them. Same seed, so the same
# species held out and the same training database. pe and se only.
set -uo pipefail
B=$HOME/bench071
W=$B/world
S=$HOME/protal-denoise/src/scripts
BIN=$HOME/protal-denoise/build
PY=$HOME/micromamba/envs/protal-db-build/bin/python
OUT=$HOME/denoise_bench
mkdir -p $OUT/logs
common=(--gtdb $W/release_p --protal $BIN/protal --simulator $BIN/simulate_metagenomes -t 6 --seed 1
        --extra-genomes $W/full/simulation/genomes_nonreps --no-binary-check --read-types pe,se --samples 8
        --read-pairs 1000,20000,200000,1000000:4 --test-samples 4 --test-read-pairs 500,10000,100000,1000000:2
        --evaluation basic)
for name in ${PIPELINES:-new base}; do
  extra=()
  [ $name = base ] && extra=(--congeners 0 --features normalized+adjacency --call-mode curve)
  [ -s $OUT/$name/model_logs/summary.txt ] && { echo "$name done"; continue; }
  echo "$(date +%T) pipeline $name"
  /usr/bin/time -v -o $OUT/logs/$name.time $PY $S/build_gtdb_database.py "${common[@]}" "${extra[@]}" --outdir $OUT/$name \
    > $OUT/logs/$name.log 2>&1 || { echo "pipeline $name failed"; tail -30 $OUT/logs/$name.log; exit 1; }
  grep -E "Elapsed|Maximum resident" $OUT/logs/$name.time
  cat $OUT/$name/model_logs/summary.txt
done
