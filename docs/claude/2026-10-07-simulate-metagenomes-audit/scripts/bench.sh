#!/usr/bin/env bash
# Before/after benchmark of simulate_metagenomes (this report): three binaries (BIN/sim_before: c9b46f9; sim_A: the exact
# changes and the genome store; sim_B: A plus the per-pair changes) on the 400 synthetic genomes of real size, every run
# on 4 pinned cores (taskset -c 0-3, nice 5), each variant REPS times, the variants alternated within a repetition.
# Prints one line per run: case, variant, repetition, wall s, CPU s (user + system), max RSS MB.
# Usage: bench.sh BIN TABLE WORK [REPS]
set -uo pipefail
bin=$1 table=$2 work=$3 reps=${4:-3}
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$work"
store="$work/store"
pin=(taskset -c 0-3 nice -n 5)
pe=(--genome_table "$table" -n 1 --species_per_sample 400 --read_length 150 --sequencer HSXt --fragment_mean 350
    --fragment_stdev 50 --seed 3 --distribution poisson_lognormal --pln_sigma 1.5 --reads_compression zstd)

timed() {  # case variant rep command...
  local name=$1 variant=$2 rep=$3; shift 3
  /usr/bin/time -f "%e %U %S %M" -o "$work/time.txt" "${pin[@]}" "$@" > "$work/last.log" 2>&1 || {
    echo "FAILED: $name $variant"; tail -3 "$work/last.log"; }
  read -r wall user sys rss < "$work/time.txt"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$variant" "$rep" "$wall" "$(echo "$user + $sys" | bc)" "$((rss / 1024))"
}

piped() {  # case variant rep sim extra...: the sample's R1 and R2 into named pipes, each drained by cat
  local name=$1 variant=$2 rep=$3 sim=$4; shift 4
  local out="$work/pipe"
  rm -rf "$out" && mkdir -p "$out/reads"
  mkfifo "$out/reads/sample_1_R1.fq.zst" "$out/reads/sample_1_R2.fq.zst"
  cat "$out/reads/sample_1_R1.fq.zst" > /dev/null & local c1=$!
  cat "$out/reads/sample_1_R2.fq.zst" > /dev/null & local c2=$!
  timed "$name" "$variant" "$rep" "$sim" "${pe[@]}" -o "$out" --total_read_pairs 2000000 -t 1 "$@"
  wait $c1 $c2
}

python3 "$here/make_inputs.py" long "$work" "$table" 400000 1
mv "$work/long_samples.tsv" "$work/tiny_samples.tsv"
python3 "$here/make_inputs.py" long "$work" "$table" 250000000 1
python3 "$here/make_inputs.py" big "$work" "$table" 50000
long() {  # case variant rep sim samples setup extra...
  local name=$1 variant=$2 rep=$3 sim=$4 samples=$5 setup=$6; shift 6
  local out="$work/long_out"
  rm -rf "$out" && mkdir -p "$out"
  sed "s#@OUT@#$out#" "$samples" > "$work/run_samples.tsv"
  timed "$name" "$variant" "$rep" "$sim" --long_samples "$work/run_samples.tsv" --long_genomes "$work/long_genomes.tsv" \
      --long_setup "$setup" --long_stats "$work/long_stats.tsv" -t 1 "$@"
}

echo -e "case\tvariant\trep\twall_s\tcpu_s\tmax_rss_mb"
for rep in $(seq 1 "$reps"); do
  rm -rf "$store"
  # per genome: 4,000 pairs of 400 genomes on one thread; the store written by the first run, read by the others
  timed shallow before "$rep" "$bin/sim_before" "${pe[@]}" -o "$work/s" --total_read_pairs 4000 -t 1
  timed shallow A "$rep" "$bin/sim_A" "${pe[@]}" -o "$work/s" --total_read_pairs 4000 -t 1
  timed shallow A_store_written "$rep" "$bin/sim_A" "${pe[@]}" -o "$work/s" --total_read_pairs 4000 -t 1 --genome_store "$store"
  timed shallow A_store_read "$rep" "$bin/sim_A" "${pe[@]}" -o "$work/s" --total_read_pairs 4000 -t 1 --genome_store "$store"
  timed shallow B "$rep" "$bin/sim_B" "${pe[@]}" -o "$work/s" --total_read_pairs 4000 -t 1
  timed shallow B_store_read "$rep" "$bin/sim_B" "${pe[@]}" -o "$work/s" --total_read_pairs 4000 -t 1 --genome_store "$store"
  # per pair: 2M pairs on one thread and on four
  for v in before A B; do
    timed deep_t1 "$v" "$rep" "$bin/sim_$v" "${pe[@]}" -o "$work/d" --total_read_pairs 2000000 -t 1
  done
  timed deep_t1 B_store_read "$rep" "$bin/sim_B" "${pe[@]}" -o "$work/d" --total_read_pairs 2000000 -t 1 --genome_store "$store"
  for v in before B; do
    timed deep_t4 "$v" "$rep" "$bin/sim_$v" "${pe[@]}" -o "$work/d" --total_read_pairs 2000000 -t 4
    timed se_t1 "$v" "$rep" "$bin/sim_$v" "${pe[@]}" -o "$work/d" --total_read_pairs 2000000 -t 1 --first_reads_only
  done
  # into named pipes: compressed as the names say, or plain
  piped pipes_t1 B_compressed "$rep" "$bin/sim_B"
  piped pipes_t1 B_plain "$rep" "$bin/sim_B" --plain_pipes
  # long reads: genome reads (a 1 kb read per genome or so) and a HiFi sample of 250 Mb
  for v in before A; do
    long hifi_genomes "$v" "$rep" "$bin/sim_$v" "$work/tiny_samples.tsv" hifi:1000:0:3
  done
  long hifi_genomes A_store_read "$rep" "$bin/sim_A" "$work/tiny_samples.tsv" hifi:1000:0:3 --genome_store "$store"
  for v in before A; do
    long hifi_250mb "$v" "$rep" "$bin/sim_$v" "$work/long_samples.tsv" hifi:15000:3000:3
  done
  long hifi_250mb A_store_read "$rep" "$bin/sim_A" "$work/long_samples.tsv" hifi:15000:3000:3 --genome_store "$store"
  # the design: 20 samples of 100-300 species from 50,000 genomes, no reads
  for v in before A; do
    timed design_50k "$v" "$rep" "$bin/sim_$v" --genome_table "$work/big_table.tsv" -o "$work/design" -n 20 \
        --species_per_sample 100-300 --strains_per_species 0.5,0.2 --congener_groups 0.25:2-5 --seed 9 --test
  done
done
du -sh "$store" | sed 's/^/store: /'
