#!/bin/bash
# The small real world's database, training and test tables and models (build_gtdb_database.py on make_inputs.sh's
# folder), scaled to 4 cores: paired-end reads (one setup, five depths, 60 samples each) and HiFi reads (three depths,
# 24 samples each), 15-60 species per sample with half of them in groups of 2-5 congeners and 2 archaeal species, 30%
# of the species and 3 genera held out of the training database (drawn by --seed), an independent test set of 15
# samples per point; no scenarios, no error reads. The column weights from GTDB's alignments and tree (the release has
# them). protal, the simulator and the scripts of WSL ~/bidx (6f8e6bb). usage: run_build.sh <inputs> <outdir> <seed>
set -eu
IN=$1; OUT=$2; SEED=$3
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
taskset -c 0-3 nice -n 5 python3 ~/bidx/src/scripts/build_gtdb_database.py --inputs $IN --outdir $OUT \
    --protal ~/bidx/build/protal --simulator ~/bidx/build/simulate_metagenomes -t 4 --seed $SEED \
    --samples 60 --read-pairs 2000,10000,50000,200000,1000000 --read-setups 150:HSXt:350:50 \
    --species-per-sample 15-60 --archaea 2 --congeners 0.5:2-5 \
    --read-types pe,pb --long-read-bases 1000000,10000000,100000000 --long-read-samples 24 \
    --test-samples 15 --test-read-pairs 10000,100000,1000000 --test-species-per-sample 15-60 \
    --test-long-read-bases 3000000,30000000 --test-long-read-samples 12 \
    --scenarios none --holdout-clades genus:3 --holdout-max-share 0.1 \
    --error-reads none --evaluation basic --no-binary-check --release 226
