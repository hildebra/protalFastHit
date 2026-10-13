#!/bin/bash
# The reduced database (build_gtdb_database.py on make_inputs.sh's folder), scaled to 4 cores as the small real
# world's (docs/claude/2026-10-10-real-small-world/run_build.sh), with Nanopore reads too: paired-end reads (one setup,
# five depths, 40 samples each), HiFi and Nanopore reads (three depths, 16 samples each), 15-60 species per sample with
# half of them in groups of 2-5 congeners and 2 archaeal species, 30% of the species and 3 genera held out of the
# training database, an independent test set of 10 samples per point; the database's default sparse-first gates, no
# models at other gates, no scenarios, no error reads.
# usage: run_build.sh <protal clone> <inputs> <outdir>     (WSL; 4 cores)
set -eu
SRC=$1; IN=$2; OUT=$3
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
    --error-reads none --evaluation basic --release 226
