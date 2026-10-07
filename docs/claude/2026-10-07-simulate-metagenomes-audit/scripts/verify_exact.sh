#!/usr/bin/env bash
# Whether two simulate_metagenomes binaries write the same bytes: paired-end samples (zstd and BGZF, with strains and
# a range of species), first reads only, long reads (HiFi, Ultima) and designs from a 50,000-genome table (--test); the
# second binary also with a genome store, when it writes the store and when it reads it.
# Usage: verify_exact.sh SIM_A SIM_B TABLE WORK   (TABLE: name, lineage, FASTA, length)
set -uo pipefail
a=$1 b=$2 table=$3 work=$4
here=$(cd "$(dirname "$0")" && pwd)
rm -rf "$work" && mkdir -p "$work"
store="$work/store"
fails=0
same() {  # label dir1 dir2: every file of dir1 the same in dir2 (but run_params.tsv, which names the command; tables
          # compared with each run's folder in their paths replaced)
  local diffs
  diffs=$(cd "$2" && find . -type f ! -name run_params.tsv | sort | while read -r f; do
    case "$f" in
      *.tsv|*.meta) cmp -s <(sed "s#$2#DIR#g" "$2/$f") <(sed "s#$3#DIR#g" "$3/$f") || echo "$f" ;;
      *) cmp -s "$2/$f" "$3/$f" || echo "$f" ;;
    esac
  done)
  if [ -z "$diffs" ] && [ "$(cd "$2" && find . -type f | wc -l)" -eq "$(cd "$3" && find . -type f | wc -l)" ]; then
    echo "same: $1 ($(cd "$2" && find . -type f | wc -l) files)"
  else
    echo "DIFFERENT: $1: ${diffs:-file count}"; fails=$((fails + 1))
  fi
}
pe() {  # sim out extra...
  local sim=$1 out=$2; shift 2
  taskset -c 0-3 nice -n 5 "$sim" --genome_table "$table" -o "$out" -n 3 --species_per_sample 30-60 \
      --strains_per_species 0.5,0.2 --total_read_pairs 30000,50000 --read_length 150 --sequencer HSXt \
      --fragment_mean 350 --fragment_stdev 50 --seed 5 --pln_sigma 1.3,2.0 -t 4 "$@" > "$out.log" 2>&1 || echo "FAILED: $out"
}
pe "$a" "$work/pe_a" --reads_compression zstd
pe "$b" "$work/pe_b" --reads_compression zstd
same "pe zstd" "$work/pe_a" "$work/pe_b"
pe "$a" "$work/gz_a" --reads_compression bgzf
pe "$b" "$work/gz_b" --reads_compression bgzf
same "pe bgzf" "$work/gz_a" "$work/gz_b"
pe "$b" "$work/store1_b" --reads_compression zstd --genome_store "$store"
same "pe, the store written" "$work/pe_a" "$work/store1_b"
pe "$b" "$work/store2_b" --reads_compression zstd --genome_store "$store"
same "pe, the store read" "$work/pe_a" "$work/store2_b"
echo "store: $(ls "$store" | wc -l) files, $(du -sh "$store" | cut -f1)"
pe "$a" "$work/se_a" --reads_compression zstd --first_reads_only
pe "$b" "$work/se_b" --reads_compression zstd --first_reads_only
same "first reads only" "$work/se_a" "$work/se_b"

python3 "$here/make_inputs.py" long "$work" "$table" 20000000 2
for setup in hifi:15000:3000:3 ultima:300:40:25:2; do
  kind=${setup%%:*}
  for run in a b b_store; do
    sim=$a; [ "$run" != a ] && sim=$b
    out="$work/long_${kind}_$run"; mkdir -p "$out"
    sed "s#@OUT@#$out#" "$work/long_samples.tsv" > "$out.samples.tsv"
    extra=(); [ "$run" = b_store ] && extra=(--genome_store "$store")
    taskset -c 0-3 nice -n 5 "$sim" --long_samples "$out.samples.tsv" --long_genomes "$work/long_genomes.tsv" \
        --long_setup "$setup" --long_stats "$out.stats.tsv" -t 4 "${extra[@]}" > "$out.log" 2>&1 || echo "FAILED: $out"
  done
  same "long $kind" "$work/long_${kind}_a" "$work/long_${kind}_b"
  same "long $kind from the store" "$work/long_${kind}_a" "$work/long_${kind}_b_store"
  cmp -s "$work/long_${kind}_a.stats.tsv" "$work/long_${kind}_b.stats.tsv" && echo "same: long $kind counts" ||
      { echo "DIFFERENT: long $kind counts"; fails=$((fails + 1)); }
done

python3 "$here/make_inputs.py" big "$work" "$table" 50000
for run in a b; do
  sim=$a; [ "$run" = b ] && sim=$b
  /usr/bin/time -f "%e s wall, %U s user, %M kB" taskset -c 0-3 nice -n 5 "$sim" --genome_table "$work/big_table.tsv" \
      -o "$work/design_$run" -n 20 --species_per_sample 100-300 --strains_per_species 0.5,0.2 \
      --congener_groups 0.25:2-5 --taxon d__Bacteria:5 --pick_random_demand_if_fail --seed 9 --test > "$work/design_$run.log" 2>&1 \
      || echo "FAILED: design $run"
  tail -1 "$work/design_$run.log"
done
same "designs (--test, 50,000 genomes, 20 samples)" "$work/design_a" "$work/design_b"
echo "$fails differences"
