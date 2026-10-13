#!/usr/bin/env bash
# What --flex_indels costs in instructions (docs/claude/2026-10-13-flex-indels): callgrind on the first READS reads of an
# ONT sample, one thread, without and with the option; the run's instructions and those of the seeding's lookups
# (ChainAnchorFinder::FindSeeds, inclusive: the lookups with the scan inlined, the seeds), and the lookups' flex cells
# (the seeding line). Writes WORK/speed.tsv.
#   bash fi_speed.sh PROTAL DATABASE ONT_READS WORK
set -uo pipefail
protal=$1 db=$2 reads=$3 work=$4
n=${READS:-2000}
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
zstdcat "$reads" | head -n $((4 * n)) > "$work/reads.fq"
echo -e "arm\tinstructions\tseeding_instructions\tcells\tseeding_instructions_per_cell" > "$work/speed.tsv"
for arm in plain flex; do
  extra=(); [ $arm = flex ] && extra=(--flex_indels ont)
  rm -rf "$work/run_$arm" "$work/callgrind.$arm"
  "${pin[@]}" valgrind --tool=callgrind --callgrind-out-file="$work/callgrind.$arm" "$protal" --db "$db" -1 "$work/reads.fq" \
      --read_type ont --prefix s -o "$work/run_$arm" -t 1 --no_strains --no_profile --sparse_first 0 "${extra[@]}" > "$work/$arm.log" 2>&1
  total=$(callgrind_annotate "$work/callgrind.$arm" 2> /dev/null | awk '/PROGRAM TOTALS/ { gsub(/,/, "", $1); print $1; exit }')
  lookup=$(callgrind_annotate --inclusive=yes "$work/callgrind.$arm" 2> /dev/null |
           awk '/ChainAnchorFinder<.*>::FindSeeds\(/ { v = $1; gsub(/,/, "", v); s += v } END { print s + 0 }')
  cells=$(grep -m1 " seeding: " "$work/$arm.log" | sed -E 's/.* blocks scanned with ([0-9]+) flex cells.*/\1/')
  echo -e "$arm\t$total\t$lookup\t$cells\t$(awk -v l="$lookup" -v c="$cells" 'BEGIN { printf "%.2f", c ? l / c : 0 }')" >> "$work/speed.tsv"
done
column -t "$work/speed.tsv"
