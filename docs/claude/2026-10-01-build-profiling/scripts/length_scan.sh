#!/bin/bash
# What simulate_metagenomes spends before its first ART call: the genome table's genomes read for their
# lengths (build_length_cache), every run. --test designs the samples without reads, so its time is that
# scan (and the design); once with the table as build_gtdb_database.py writes it (no lengths), once with a
# genome_length column, which the simulator takes instead of reading the genomes.
# usage: length_scan.sh GENOME_TABLE OUT
set -euo pipefail
table=$1 out=$2
sim=${SIM:-$HOME/build-prof/build/simulate_metagenomes}
mkdir -p "$out"
# The table with lengths: the padded copies' sizes from pad_genomes.sh's plan.
plan=$(dirname "$table")/plan.tsv
awk -F'\t' -v OFS='\t' 'NR == FNR { len[$1] = $3; next } { print $1, $2, $3, len[$1] }' "$plan" "$table" > "$out/with_lengths.tsv"
for variant in plain with_lengths; do
  t=$([ $variant = plain ] && echo "$table" || echo "$out/with_lengths.tsv")
  /usr/bin/time -f "$variant: wall %e s, user %U s, sys %S s" -a -o "$out/times.txt" \
    "$sim" --genome_table "$t" -o "$out/$variant" -n 2 --total_read_pairs 1000 --species_per_sample 20-200 \
    --seed 1 --test > "$out/$variant.log" 2>&1
done
cat "$out/times.txt"
