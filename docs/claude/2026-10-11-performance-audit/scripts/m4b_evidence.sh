#!/usr/bin/env bash
# The cohort's kept memory with and without each sample's record evidence dropped after its profile is written
# (tim: PROTAL_AUDIT_DROP_EVIDENCE), 1 and 24 samples, strains on, no qcMSA, 4 threads on cores 0-3.
set -uo pipefail
W=$HOME/perf9; S=$W/samples/cohort/reads; R=$W/runs/cohort_evidence; DB=$W/arc/protal_db
mkdir -p "$R/logs"
printf 'variant\tsamples\twall_s\tmax_rss_mb\tafter_aligning\tafter_profiling\n' > "$R/evidence.tsv"
mem() { grep -m1 "^Memory after $1:" "$2" | sed -E 's/.*: ([0-9.]+) GB resident, peak ([0-9.]+) GB/\1\/\2/'; }
for rep in 1 2; do
  for n in 1 24; do
    for variant in kept dropped; do
      out=$R/work; rm -rf "$out"; mkdir -p "$out"; map=$out/samples.map
      { printf '#OUTPUT_DIR\t%s\n' "$out"; printf '#SAMPLEID\tFIRST\tSECOND\tPREFIX\tREAD_TYPE\n'
        for i in $(seq 1 "$n"); do printf 'c%s\t%s\t%s\tc%s\tpe\n' "$i" "$S/c_${i}_R1.fq.gz" "$S/c_${i}_R2.fq.gz" "$i"; done; } > "$map"
      env=(); [ "$variant" = dropped ] && env=(PROTAL_AUDIT_DROP_EVIDENCE=1)
      log=$R/logs/$variant.$n.$rep.log
      /usr/bin/time -f "%e %M" -o "$out/time.txt" env "${env[@]}" taskset -c 0-3 nice -n 5 "$W/tim/build/protal" --db "$DB" --map "$map" \
          -t 4 --no_qcmsa > "$log" 2>&1 || echo "FAILED $variant $n"
      read -r wall rss < <(tail -1 "$out/time.txt")
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$variant" "$n" "$wall" "$((rss / 1024))" "$(mem aligning "$log")" "$(mem profiling "$log")" >> "$R/evidence.tsv"
      # the profiles must not change: keep their checksums
      (cd "$out/profiles" 2> /dev/null && cat *.profile | md5sum | sed "s/-/$variant.$n.$rep/") >> "$R/profiles.md5"
      rm -rf "$out"
    done
  done
done
echo DONE
