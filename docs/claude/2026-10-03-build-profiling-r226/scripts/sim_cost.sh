#!/bin/bash
# What a paired-end sample of simulate_metagenomes costs, per genome and per read: one simulator run of
# SAMPLES samples at each of READ_PAIRS (one thread, the build's species design, 150 bp HSXt), every ART call
# timed (timed_tool.sh of the 2026-10-01 report). Prints, per run: wall, CPU, genomes simulated, ART calls and
# their wall/CPU, and what is left (the simulator: genome decompression, the reads' BGZF compression, bookkeeping).
#
# usage: sim_cost.sh OUT
# needs: ~/bprof/b2/db/genomes.tsv (the padded real-sized genomes with lengths), SIM (simulate_metagenomes)
set -euo pipefail
out=$1
here=$(cd "$(dirname "$0")" && pwd)
old=$(cd "$here/../../2026-10-01-build-profiling/scripts" && pwd)
SIM=${SIM:-$HOME/protal-fp/build/simulate_metagenomes}
TABLE=${TABLE:-$HOME/bprof/b2/db/genomes.tsv}
SAMPLES=${SAMPLES:-2}
READ_PAIRS=${READ_PAIRS:-1000 500000}
mkdir -p "$out"
bash "$old/make_wrappers.sh" "$out/bin" "$out/tool_times.tsv"
export PATH="$out/bin:$PATH"
printf 'pairs\tsamples\twall_s\tcpu_s\tgenomes\tart_calls\tart_wall_s\tart_cpu_s\trest_wall_s\tmax_rss_kb\n' > "$out/summary.tsv"
for pairs in $READ_PAIRS; do
  run=$out/p$pairs
  rm -rf "$run"; : > "$out/tool_times.tsv"
  /usr/bin/time -f '%e\t%U\t%S\t%M' -o "$out/p$pairs.time" "$SIM" --genome_table "$TABLE" -o "$run" -n "$SAMPLES" \
    --sample_prefix "p${pairs}_s" --total_read_pairs "$pairs" --species_per_sample 20-200 --read_length 150 \
    --sequencer HSXt --fragment_mean 350 --fragment_stdev 50 --seed 7 -t 1 --strains_per_species 0.3,0.1 \
    --taxon d__Archaea:2 --pick_random_demand_if_fail > "$out/p$pairs.log" 2>&1 \
    || { echo "simulator failed, see $out/p$pairs.log"; exit 1; }
  genomes=$(tail -n +2 "$run/manifest.tsv" | wc -l)
  awk -v pairs="$pairs" -v samples="$SAMPLES" -v genomes="$genomes" -v t="$(cat "$out/p$pairs.time")" '
    BEGIN { split(t, a, "\t") }
    $1 == "art_illumina" { n++; w += $2; c += $3 + $4 }
    END { printf "%s\t%s\t%.1f\t%.1f\t%d\t%d\t%.1f\t%.1f\t%.1f\t%s\n", pairs, samples, a[1], a[2] + a[3], genomes, n, w, c, a[1] - w, a[4] }' \
    "$out/tool_times.tsv" >> "$out/summary.tsv"
  cp "$out/tool_times.tsv" "$out/p$pairs.tool_times.tsv"
done
# One genome by itself: decompressing it, ART at a few pairs and at many.
genome=$(awk -F'\t' 'NR == 1 { print $3 }' "$TABLE")
{
  printf 'step\twall_s\n'
  s=$(date +%s.%N); zcat "$genome" > "$out/one.fna"; e=$(date +%s.%N); printf 'decompress one genome (%s bytes)\t%.3f\n' "$(stat -c %s "$out/one.fna")" "$(echo "$e - $s" | bc)"
  for cov in 0.001 1; do
    s=$(date +%s.%N); art_illumina -ss HSXt -i "$out/one.fna" -p -l 150 -f $cov -m 350 -s 50 -na -rs 1 -o "$out/one_" > /dev/null 2>&1; e=$(date +%s.%N)
    printf 'ART at coverage %s (%s pairs)\t%.3f\n' "$cov" "$(( $(wc -l < "$out/one_1.fq") / 4 ))" "$(echo "$e - $s" | bc)"
  done
} > "$out/one_genome.tsv"
rm -f "$out/one.fna" "$out"/one_*.fq
cat "$out/summary.tsv" "$out/one_genome.tsv"
