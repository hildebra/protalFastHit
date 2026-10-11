#!/usr/bin/env bash
# After impl_asan.sh: the timed runs and the profile-only runs once more, with the final build (binary headers read
# without their frames): runs3/runs.tsv, cohort.tsv and tables3/profile_only.tsv. 4 pinned cores, REPS rounds.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
until grep -q "ASAN DONE" "$A/impl_asan.out" 2> /dev/null; do sleep 15; done
W=$HOME/impl/perf; reps=${REPS:-3}
REF=$HOME/impl/ref/build/protal; WORK=$HOME/impl/work/build/protal
echo "[$(date +%T)] timed runs (runs3)"
OUT=$W/runs3 bash "$A/impl_perf_runs.sh" 2>&1 | tail -30
echo "[$(date +%T)] profile-only runs (tables3)"
R=$W/tables3; rm -rf "$R"; mkdir -p "$R/logs"
pin=(taskset -c 0-3 nice -n 5)
secs() { awk '{ t = 0; for (i = 1; i <= NF; i++) { if ($i ~ /^[0-9]+m$/) t += 60 * $i; else if ($i ~ /^[0-9]+s$/) t += $i; else if ($i ~ /^[0-9]+ms$/) t += $i / 1000 } printf "%.3f", t }'; }
mem() { grep -m1 "^Memory after $1:" "$2" | sed -E 's/.*: ([0-9.]+) GB resident, peak ([0-9.]+) GB/\1\/\2/'; }
sam=$(find "$W/post/ref_full" -name "*.sam.zst" | head -1)
printf 'variant\trep\twall_s\tmax_rss_mb\ttables_s\tpreload_s\tprofiling_s\tafter_profiling\n' > "$R/profile_only.tsv"
for rep in $(seq 1 "$reps"); do
  for v in ref_text work_text work_bin; do
    case $v in ref_text) bin=$REF db=$W/db_ref ;; work_text) bin=$WORK db=$W/db_ref ;; work_bin) bin=$WORK db=$W/db_bin ;; esac
    out=$R/out; rm -rf "$out"
    log=$R/logs/po.$v.$rep.log
    /usr/bin/time -f "%e %M" -o "$R/time.txt" "${pin[@]}" "$bin" --db "$db/database.protal" --profile_only "$sam" -o "$out" -t 4 --no_qcmsa \
        > "$log" 2>&1 || echo "FAILED $v $rep"
    read -r wall rss < <(tail -1 "$R/time.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$v" "$rep" "$wall" "$((rss / 1024))" \
      "$(grep -m1 '^Loading the taxonomy, models and tables took' "$log" | secs)" "$(grep -m1 '^Preload genomes took' "$log" | secs)" \
      "$(grep -m1 '^Profiling took' "$log" | secs)" "$(mem profiling "$log")" >> "$R/profile_only.tsv"
    rm -rf "$out"
  done
done
cat "$R/profile_only.tsv"
echo "[$(date +%T)] RUNS3 DONE"
