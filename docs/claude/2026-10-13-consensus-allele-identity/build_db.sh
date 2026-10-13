#!/bin/bash
# The reduced database of the mock local test (../2026-10-13-mock-local-test/run_build.sh: 761 r226 species, 4 cores,
# seed 1), built by one protal clone, with more build_gtdb_database.py options (e.g. --features).
# usage: build_db.sh <protal clone> <inputs> <outdir> [build_gtdb_database.py options ...]     (WSL; 4 cores)
set -eu
SRC=$1; IN=$2; OUT=$3; shift 3
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
taskset -c 0-3 nice -n 5 python3 $SRC/scripts/build_gtdb_database.py --inputs $IN --outdir $OUT \
    --protal $SRC/build/protal --simulator $SRC/build/simulate_metagenomes -t 4 --seed 1 \
    --samples 40 --read-pairs 2000,10000,50000,200000,1000000 --read-setups 150:HSXt:350:50 \
    --species-per-sample 15-60 --archaea 2 --congeners 0.5:2-5 \
    --read-types pe,pb,ont --long-read-bases 1000000,10000000,100000000 --long-read-samples 16 \
    --pbsim-models $HOME/micromamba/envs/protal-db-build/data \
    --test-samples 10 --test-read-pairs 10000,100000,1000000 --test-species-per-sample 15-60 \
    --test-long-read-bases 3000000,30000000 --test-long-read-samples 8 \
    --scenarios none --holdout-clades genus:3 --holdout-max-share 0.1 --model-gates none \
    --error-reads none --evaluation basic --release 226 "$@"
