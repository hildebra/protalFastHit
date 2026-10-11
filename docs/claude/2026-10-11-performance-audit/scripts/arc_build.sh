#!/usr/bin/env bash
# The audit's local world: every GTDB r226 archaeal species (6,968; markers of all 17,245 genomes) and the 3,259 whole
# genomes of ~/archaea_inputs, built and trained by 0.8.1's own pipeline (build_gtdb_database.py at b3ac271) with a
# small design, on 4 pinned cores. Its protal_db is a production database at 1/36 of r226's genes.
set -uo pipefail
W=$HOME/perf9
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
start=$(date +%s)
taskset -c 0-3 nice -n 5 python3 "$W/v081/src/scripts/build_gtdb_database.py" --gtdb "$HOME/archaea_inputs/release" --extra-genomes "$HOME/archaea_inputs/genomes" --simulate-species "$HOME/archaea_inputs/simulation_species.txt" --outdir "$W/arc" \
    --protal "$W/v081/build/protal" --simulator "$W/v081/build/simulate_metagenomes" -t 4 --seed 1 \
    --samples 24 --read-pairs 10000,200000,1000000 --read-setups 150:HSXt:350:50 --species-per-sample 30-150 \
    --read-types pe,pb --long-read-bases 10000000,100000000 --long-read-samples 12 \
    --test-samples 6 --test-read-pairs 100000,1000000 --test-species-per-sample 30-150 \
    --test-long-read-bases 30000000 --test-long-read-samples 6 --scenarios none --holdout-clades genus:3 \
    --holdout-max-share 0.1 --error-reads none --evaluation basic --no-binary-check --release 226 2>&1
echo "EXIT $? after $(( $(date +%s) - start )) s"
