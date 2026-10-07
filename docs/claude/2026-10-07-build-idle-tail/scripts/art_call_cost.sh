#!/usr/bin/env bash
# The cost of one simulate_metagenomes genome step, on one core: the genome inflated to a plain copy, and one
# art_illumina call (HSXt, 150 bp pairs, fragments 350/50, no alignment files) at a soil_shallow genome's few pairs
# and at many pairs (the per-pair cost).
# Usage: art_call_cost.sh GENOME_DIR WORK_DIR [ART]   (GENOME_DIR: *.fna.gz; the first 10 are used)
set -euo pipefail
genomes=$1 work=$2 art=${3:-art_illumina}
mkdir -p "$work"
ls "$genomes"/*.fna.gz | head -10 > "$work/genomes.txt"
t() { local s=$(date +%s.%N); "$@" > /dev/null 2>&1; echo "$(date +%s.%N) - $s" | bc; }
total() { awk '{s += $1} END {printf "%.1f ms per genome\n", s / NR * 1000}'; }
echo "inflate (gzip -dc to a file):"
i=0; while read -r g; do i=$((i + 1)); t sh -c "gzip -dc '$g' > '$work/g$i.fna'"; done < "$work/genomes.txt" | total
len() { grep -v '>' "$1" | tr -d '\n' | wc -c; }
for pairs in 100 400 20000; do
  echo "art_illumina at $pairs pairs:"
  i=0; while read -r g; do i=$((i + 1)); L=$(len "$work/g$i.fna"); f=$(echo "$pairs * 300 / $L" | bc -l)
    t "$art" -ss HSXt -i "$work/g$i.fna" -p -l 150 -f "$f" -m 350 -s 50 -na -rs $i -o "$work/r$i"; done < "$work/genomes.txt" | total
done
echo "art_illumina's start alone (a 1 kb reference, ~3 pairs):"
head -c 1100 "$work/g1.fna" | head -n 15 > "$work/tiny.fna"
for i in $(seq 10); do t "$art" -ss HSXt -i "$work/tiny.fna" -p -l 150 -f 1 -m 350 -s 50 -na -rs $i -o "$work/tiny"; done | total
rm -f "$work"/g*.fna "$work"/r*.fq "$work"/tiny*
