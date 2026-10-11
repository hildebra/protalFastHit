#!/usr/bin/env bash
# The audit's short-read samples from the archaea world's genome table (real genomes): one 1M-pair community (its R1 is
# the single-end sample) and a cohort of 24 samples of 50,000 pairs. 4 pinned cores.
set -uo pipefail
S=$HOME/perf9/samples
SIM=$HOME/perf9/v081/build/simulate_metagenomes
T=$HOME/perf9/arc/work/genomes.tsv
mkdir -p "$S"
pin=(taskset -c 0-3 nice -n 5)
"${pin[@]}" "$SIM" --genome_table "$T" -o "$S/pe1m" -n 1 --sample_prefix pe1m --total_read_pairs 1000000 \
    --species_per_sample 120 --strains_per_species 0.3 --congener_groups 0.25:2-5 --distribution lognormal --pln_sigma 1.5 \
    --read_length 150 --sequencer HSXt --fragment_mean 350 --fragment_stdev 50 --seed 11 -t 4 --reads_compression bgzf \
    > "$S/pe1m.log" 2>&1 || { echo "FAIL pe1m"; tail -5 "$S/pe1m.log"; }
"${pin[@]}" "$SIM" --genome_table "$T" -o "$S/cohort" -n 24 --sample_prefix c --total_read_pairs 50000 \
    --species_per_sample 30-150 --strains_per_species 0.3 --congener_groups 0.25:2-5 --distribution lognormal --pln_sigma 1.5 \
    --read_length 150 --sequencer HSXt --fragment_mean 350 --fragment_stdev 50 --seed 21 -t 4 --reads_compression bgzf \
    > "$S/cohort.log" 2>&1 || { echo "FAIL cohort"; tail -5 "$S/cohort.log"; }
find "$S" -maxdepth 3 | head -40
du -sh "$S"/*
find ~/micromamba -name 'QSHMM-ONT-HQ.model' 2> /dev/null | head -2
