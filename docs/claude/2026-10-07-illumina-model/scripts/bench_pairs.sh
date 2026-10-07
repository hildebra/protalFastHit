#!/usr/bin/env bash
# simulate_metagenomes's paired-end reads (IlluminaSimulator) on the 400 synthetic genomes of real size: CPU per read
# pair (a deep sample) and per genome of a sample (a shallow one: mostly reading the genomes), at 1 and 4 threads.
# Usage: bench_pairs.sh SIMULATOR GENOME_DIR WORK_DIR
set -euo pipefail
sim=$1 genomes=$2 work=$3
mkdir -p "$work"
table="$work/genomes.tsv"
if [ ! -s "$table" ]; then  # name, taxonomy, FASTA, length (the simulator reads no genome for its length)
  i=0
  for f in "$genomes"/*.fna.gz; do
    i=$((i + 1))
    len=$(zcat "$f" | grep -v '>' | tr -d '\n' | wc -c)
    printf 'G%d\td__Bacteria;p__P;c__C;o__O;f__F;g__G%d;s__G%d sp\t%s\t%s\n' $i $i $i "$f" $len >> "$table"
  done
fi
run() {  # label pairs species threads
  local out="$work/$1_t$4"
  rm -rf "$out"
  /usr/bin/time -f "%e %U %S %M" -o "$work/time.txt" "$sim" --genome_table "$table" -o "$out" -n 1 --species_per_sample "$3" \
      --total_read_pairs "$2" --read_length 150 --sequencer HSXt --fragment_mean 350 --fragment_stdev 50 --seed 3 \
      --distribution poisson_lognormal --pln_sigma 1.5 -t "$4" --reads_compression zstd > "$work/$1_t$4.log" 2>&1
  read -r wall user sys rss < "$work/time.txt"
  local genomes_with_pairs
  genomes_with_pairs=$(awk -F'\t' 'NR > 1 && $6 > 0' "$out/manifest.tsv" | wc -l)
  echo -e "$1\t$2\t$genomes_with_pairs\t$4\t$wall\t$(echo "$user + $sys" | bc)\t$((rss / 1024))"
}
echo -e "run\tpairs\tgenomes\tthreads\twall_s\tcpu_s\tmax_rss_mb"
run deep 2000000 400 1
run deep 2000000 400 4
run shallow 4000 400 1
