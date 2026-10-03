#!/bin/bash
# simulate_metagenomes with a sample's genomes on threads: the same reads as one thread, and as the binary before
# the change. Runs OLD -t 1, NEW -t 1, NEW -t 6 (one sample: 6 genome threads), NEW -t 4 with 2 samples (2 genome
# threads each), and NEW --test (the design run of the collector), and compares reads and manifests.
#
# usage: sim_identity.sh OUT   (OLD, NEW: simulate_metagenomes binaries; TABLE: genome table with lengths)
set -euo pipefail
out=$1
OLD=${OLD:-$HOME/protal-fp/build/simulate_metagenomes}
NEW=${NEW:-$HOME/bt26/src/build/simulate_metagenomes}
TABLE=${TABLE:-$HOME/bprof/b2/db/genomes.tsv}
PAIRS=${PAIRS:-20000}
rm -rf "$out"; mkdir -p "$out"
sim() {  # name binary threads samples [extra]
  local name=$1 bin=$2 t=$3 n=$4; shift 4
  /usr/bin/time -f "$name %e s wall %U s user" "$bin" --genome_table "$TABLE" -o "$out/$name" -n "$n" --sample_prefix x_s \
    --total_read_pairs "$PAIRS" --species_per_sample 20-200 --read_length 150 --sequencer HSXt --fragment_mean 350 \
    --fragment_stdev 50 --seed 11 -t "$t" --strains_per_species 0.3,0.1 --taxon d__Archaea:2 \
    --pick_random_demand_if_fail --protal_metafile "$out/$name/protal" "$@" > "$out/$name.log" 2> "$out/$name.err"
  tail -1 "$out/$name.err"
}
sim old_t1 "$OLD" 1 2
sim new_t1 "$NEW" 1 2
sim new_t6 "$NEW" 6 1
sim new_t4 "$NEW" 4 2
sim design "$NEW" 1 2 --test
for run in new_t1 new_t4; do
  for f in "$out"/old_t1/reads/*.fq.gz; do
    cmp -s "$f" "$out/$run/reads/$(basename "$f")" && echo "$run $(basename "$f") identical" || echo "$run $(basename "$f") DIFFERS"
  done
done
for f in "$out"/old_t1/reads/x_s_1_R*.fq.gz; do
  cmp -s "$f" "$out/new_t6/reads/$(basename "$f")" && echo "new_t6 $(basename "$f") identical" || echo "new_t6 $(basename "$f") DIFFERS"
done
strip() { cut -f1-8,11- "$1"; }  # without fastq_r1, fastq_r2
for run in new_t1 new_t4 design; do
  cmp -s <(strip "$out/old_t1/manifest.tsv") <(strip "$out/$run/manifest.tsv") && echo "$run manifest identical" || echo "$run manifest DIFFERS"
done
ls "$out"/*/ -d | while read d; do [ -z "$(ls -d "$d"/*_tmp 2>/dev/null)" ] || echo "$d left temporary folders"; done
