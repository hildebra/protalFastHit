#!/usr/bin/env bash
# What a cohort keeps per sample after profiling: the 24 cohort samples aligned once (SAMs kept), then a profile-only run
# of all of them under massif (no index loaded, so the samples' slots stand out), --no_strains --no_qcmsa, 4 threads on
# cores 0-3. ms_print's snapshots: the heap over the run and its allocation trees.
set -uo pipefail
W=$HOME/perf9; S=$W/samples/cohort/reads; R=$W/runs/massif; DB=$W/arc/protal_db
mkdir -p "$R/sams"
map=$R/align.map
{ printf '#OUTPUT_DIR\t%s\n' "$R/sams"; printf '#SAMPLEID\tFIRST\tSECOND\tPREFIX\tREAD_TYPE\n'
  for i in $(seq 1 24); do printf 'c%s\t%s\t%s\tc%s\tpe\n' "$i" "$S/c_${i}_R1.fq.gz" "$S/c_${i}_R2.fq.gz" "$i"; done; } > "$map"
ls "$R/sams/"*.sam.zst > /dev/null 2>&1 || taskset -c 0-3 nice -n 5 "$W/v081/build/protal" --db "$DB" --map "$map" -t 4 --no_profile > "$R/align.log" 2>&1
ls "$R/sams/"*.sam.zst | wc -l
rm -rf "$R/out"; mkdir -p "$R/out"
taskset -c 0-3 nice -n 5 valgrind --tool=massif --threshold=0.2 --detailed-freq=1 --max-snapshots=120 --massif-out-file="$R/massif.out" \
    "$W/v081/build/protal" --db "$DB" --profile_only "$R/sams/*.sam.zst" --read_type pe -o "$R/out" -t 4 --no_strains --no_qcmsa > "$R/massif.log" 2>&1
echo "exit $?"
ms_print --threshold=0.2 "$R/massif.out" > "$R/massif.txt"
grep -E "^Memory after|took" "$R/massif.log" | tail -8
head -60 "$R/massif.txt"
