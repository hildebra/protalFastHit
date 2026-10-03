#!/bin/bash
# A collection's simulations (collect_training_data.py --simulate_only) before and after the changes, on the padded
# real-sized genomes: pe, pb, ont; 1,000 / 20,000 / 500,000 (2 samples) / 2,000,000 (1 sample) read pairs of 150 bp,
# 3 samples per point; long reads of 0.3 / 6 / 30 / 600 (1 sample) Mb. Then compares the samples' reads (decompressed:
# Python's gzip writes the time into the header): the paired-end reads and the unchunked long-read samples must be
# the same, the 600 Mb ones (chunked) differ.
#
# usage: collect_ab.sh OUT [THREADS]
#   OLD_SCRIPTS (a folder with the collector of HEAD), OLD_SIM; NEW_SCRIPTS, NEW_SIM
set -uo pipefail
out=$1 threads=${2:-6}
OLD_SCRIPTS=${OLD_SCRIPTS:-$HOME/bt26/old_scripts}
OLD_SIM=${OLD_SIM:-$HOME/protal-fp/build/simulate_metagenomes}
NEW_SCRIPTS=${NEW_SCRIPTS:-$HOME/bt26/src/scripts}
NEW_SIM=${NEW_SIM:-$HOME/bt26/src/build/simulate_metagenomes}
TABLE=${TABLE:-$HOME/bprof/b2/db/genomes.tsv}
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
mkdir -p "$out"
run() {  # name scripts simulator
  local name=$1 scripts=$2 sim=$3
  rm -rf "$out/$name"
  uptime > "$out/$name.load"
  /usr/bin/time -f "%e s wall, %U s user, %S s system, %M kB" -o "$out/$name.time" \
    python "$scripts/collect_training_data.py" --db /nonexistent --genome_table "$TABLE" -o "$out/$name" \
    --simulator "$sim" --pbsim "$HOME/micromamba/envs/protal-db-build/bin/pbsim" --read_types pe,pb,ont \
    --read_pairs 1000,20000,500000:2,2000000:1 --read_setups 150:HSXt:350:50 --samples 3 \
    --long_read_bases 300000,6000000,30000000,600000000:1 --species_per_sample 20-200 --strains_per_species 0.3,0.1 \
    --archaea 2 --seed 1 -t "$threads" --simulate_only > "$out/$name.log" 2>&1
  echo "$name: $(cat "$out/$name.time")"
  uptime >> "$out/$name.load"
}
run before "$OLD_SCRIPTS" "$OLD_SIM"
run after "$NEW_SCRIPTS" "$NEW_SIM"
cd "$out"
for f in $(cd before && find points -name '*.fq.gz' | sort); do
  if cmp -s <(zcat "before/$f") <(zcat "after/$f"); then echo "identical $f"; else echo "DIFFERS $f"; fi
done > compare.txt
grep -c identical compare.txt; grep DIFFERS compare.txt
