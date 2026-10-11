#!/usr/bin/env bash
# A cohort's memory: 0.8.1 on maps of 1, 6 and 24 samples of 50,000 pairs (the cohort samples), strains on (the
# default) but no qcMSA, by default and with --profile_ahead (as the database builds run), and with --no_strains;
# 4 threads on 4 pinned cores. Writes cohort.tsv: variant, samples, wall, max RSS, and protal's memory lines.
set -uo pipefail
W=$HOME/perf9; S=$W/samples/cohort/reads; R=${OUT:-$W/runs/cohort}; DB=${DB:-$W/arc/protal_db}
mkdir -p "$R/logs"
printf 'variant\tsamples\twall_s\tmax_rss_mb\tafter_index\tafter_aligning\tafter_profiling\n' > "$R/cohort.tsv"
mem() { grep -m1 "^Memory after $1:" "$2" | sed -E 's/.*: ([0-9.]+) GB resident, peak ([0-9.]+) GB/\1\/\2/'; }
for n in 1 6 24; do
  for variant in default profile_ahead no_strains; do
    out=$R/work; rm -rf "$out"; mkdir -p "$out"
    map=$out/samples.map
    { printf '#OUTPUT_DIR\t%s\n' "$out"; printf '#SAMPLEID\tFIRST\tSECOND\tPREFIX\tREAD_TYPE\n'
      for i in $(seq 1 "$n"); do printf 'c%s\t%s\t%s\tc%s\tpe\n' "$i" "$S/c_${i}_R1.fq.gz" "$S/c_${i}_R2.fq.gz" "$i"; done; } > "$map"
    case $variant in default) opts=() ;; profile_ahead) opts=(--profile_ahead) ;; no_strains) opts=(--no_strains) ;; esac
    log=$R/logs/$variant.$n.log
    /usr/bin/time -f "%e %M" -o "$out/time.txt" taskset -c 0-3 nice -n 5 "$W/v081/build/protal" --db "$DB" --map "$map" -t 4 \
        --no_qcmsa "${opts[@]}" > "$log" 2>&1 || echo "FAILED $variant $n"
    read -r wall rss < "$out/time.txt"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$variant" "$n" "$wall" "$((rss / 1024))" "$(mem 'loading the index' "$log")" \
      "$(mem aligning "$log")" "$(mem profiling "$log")" >> "$R/cohort.tsv"
    rm -rf "$out"
  done
done
echo DONE
