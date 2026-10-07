#!/usr/bin/env bash
# Instructions (callgrind) of the three simulate_metagenomes binaries of bench.sh on one thread, which the laptop's
# changing core speeds do not touch: paired-end reads (100,000 pairs of 40 genomes), first reads only, a shallow sample
# (400 pairs of 40 genomes: mostly reading genomes) with and without a genome store; and the functions' shares in B.
# Usage: profile_ir.sh BIN TABLE WORK
set -uo pipefail
bin=$1 table=$2 work=$3
mkdir -p "$work"
store="$work/store"
pin=(taskset -c 0-3 nice -n 5)
pe=(--genome_table "$table" -n 1 --species_per_sample 40 --read_length 150 --sequencer HSXt --fragment_mean 350
    --fragment_stdev 50 --seed 3 --distribution poisson_lognormal --pln_sigma 1.5 --reads_compression zstd -t 1)
ir() {  # label binary args...: the program's instructions
  local label=$1 sim=$2; shift 2
  "${pin[@]}" valgrind --tool=callgrind --callgrind-out-file="$work/cg.$label" "$sim" "${pe[@]}" -o "$work/out_$label" "$@" \
      > "$work/$label.log" 2>&1 || echo "FAILED: $label"
  printf '%s\t%s\n' "$label" "$(grep '^summary:' "$work/cg.$label" | awk '{print $2}')"
}
echo -e "case\tinstructions"
for v in before A B; do ir "deep_$v" "$bin/sim_$v" --total_read_pairs 100000; done
for v in before B; do ir "se_$v" "$bin/sim_$v" --total_read_pairs 100000 --first_reads_only; done
for v in before A; do ir "shallow_$v" "$bin/sim_$v" --total_read_pairs 400; done
rm -rf "$store"
"$bin/sim_A" "${pe[@]}" -o "$work/out_fill" --total_read_pairs 400 --genome_store "$store" > /dev/null 2>&1
ir shallow_A_store_read "$bin/sim_A" --total_read_pairs 400 --genome_store "$store"
ir deep_B_store_read "$bin/sim_B" --total_read_pairs 100000 --genome_store "$store"
for label in deep_before deep_B shallow_before shallow_A; do
  echo "== $label: functions by own instructions"
  callgrind_annotate "$work/cg.$label" 2>/dev/null | grep -E '^ *[0-9,]+ \(' | head -14
done
