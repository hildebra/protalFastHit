#!/usr/bin/env bash
# The benchmark samples, by 0.7.1's training-data collector (so in the design its training uses), from all 900
# species of the world (the 135 unknown to every database included: their reads land on relatives), seeds of
# their own. Paired-end: 2x150 (HiSeq X profile) and 2x100 (HiSeq 2500), 1,000, 10,000 and 500,000 read pairs,
# 4 samples each; single-end: their first reads; PacBio and Nanopore (pbsim3): 3 and 90 Mb from the same
# communities; 10-300 species per sample, Poisson-lognormal abundances of sigma 2, a second strain in half of
# the species and a third in a fifth, 6 archaeal species per sample. Deep: 2 samples of 5M pairs 2x150. The
# collector profiles them with 0.7.1 too (against $SAMPLES_DB, by default V071's database); profile.sh runs every
# version on them.
set -uo pipefail
B=${BENCH:-$HOME/bench071}
W=$B/world
PY=${PY:-$HOME/micromamba/envs/protal-db-build/bin/python}
PBSIM=${PBSIM:-$HOME/micromamba/envs/protal-db-build/bin/pbsim}
C=$B/src/0.7.1/scripts/collect_training_data.py
common=(--db ${SAMPLES_DB:-$B/V071/protal_db} --genome_table $W/full/simulation/genomes.tsv --protal $B/bin/protal-0.7.1
        --simulator $B/bin/simulate_metagenomes --species_per_sample 10-300 --abundance lognormal:2.0
        --strains_per_species 0.5,0.2 --archaea 6 --pbsim $PBSIM -t 6)
$PY $C "${common[@]}" -o $B/samples --samples 4 --read_pairs 1000,10000,500000 \
  --read_setups 150:HSXt:350:50,100:HS20:300:40 --read_types pe,se,pb,ont --long_read_bases 3000000,90000000 \
  --seed 501 > $B/logs/samples.log 2>&1 || { echo "samples failed"; tail -20 $B/logs/samples.log; exit 1; }
$PY $C "${common[@]}" -o $B/samples_deep --samples 2 --read_pairs 5000000 --read_setups 150:HSXt:350:50 \
  --read_types pe --seed 601 > $B/logs/samples_deep.log 2>&1 || { echo "deep samples failed"; tail -20 $B/logs/samples_deep.log; exit 1; }
ls $B/samples/points $B/samples_deep/points
