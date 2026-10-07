#!/usr/bin/env bash
# The genome store's writer before and after its table-driven packing (sim_B: a switch per base; sim_B2: four bases a
# byte by a lookup table): the same store files, and the time of a shallow sample (4,000 pairs of 400 genomes, one
# thread) that writes the store, against one that reads FASTAs only. Alternated, REPS times.
# Usage: store_writer.sh BIN TABLE WORK [REPS]
set -uo pipefail
bin=$1 table=$2 work=$3 reps=${4:-3}
rm -rf "$work" && mkdir -p "$work"
pe=(--genome_table "$table" -n 1 --species_per_sample 400 --read_length 150 --sequencer HSXt --fragment_mean 350
    --fragment_stdev 50 --seed 3 --distribution poisson_lognormal --pln_sigma 1.5 --reads_compression zstd
    --total_read_pairs 4000 -t 1)
run() {  # label sim extra...
  local label=$1 sim=$2 rep=$3; shift 3
  /usr/bin/time -f "%e %U %S" -o "$work/time.txt" taskset -c 0-3 nice -n 5 "$sim" "${pe[@]}" -o "$work/out" "$@" \
      > "$work/last.log" 2>&1 || echo "FAILED: $label"
  read -r wall user sys < "$work/time.txt"
  printf '%s\t%s\t%s\t%s\n' "$label" "$rep" "$wall" "$(echo "$user + $sys" | bc)"
}
echo -e "case\trep\twall_s\tcpu_s"
for rep in $(seq 1 "$reps"); do
  run fasta_only "$bin/sim_B2" "$rep"
  rm -rf "$work/store_B" && run store_written_B "$bin/sim_B" "$rep" --genome_store "$work/store_B"
  rm -rf "$work/store_B2" && run store_written_B2 "$bin/sim_B2" "$rep" --genome_store "$work/store_B2"
done
same=0 diff=0
for f in "$work"/store_B/*.g2b; do
  if cmp -s "$f" "$work/store_B2/$(basename "$f")"; then same=$((same + 1)); else diff=$((diff + 1)); fi
done
echo "store files: $same the same, $diff different, $(ls "$work/store_B2" | wc -l) written by B2"
