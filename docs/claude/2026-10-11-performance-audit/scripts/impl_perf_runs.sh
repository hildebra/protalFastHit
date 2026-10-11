#!/usr/bin/env bash
# Before and after on the archaea world: HEAD (ref) on the database of text tables (db_ref), the changes (work) on the
# same database (the text read once, the profiling's tables beside the alignment) and on the database with the binary
# forms (db_bin); the 1M-pair pe sample and 100 Mb of HiFi, 4 threads, REPS rounds alternated; then the cohort of 24
# samples, ref against work. One line per run in runs.tsv: wall, max RSS, the start-up's timers, memory lines.
set -uo pipefail
W=$HOME/impl/perf; S=$W/samples; R=${OUT:-$W/runs}; reps=${REPS:-3}
mkdir -p "$R/logs"
pin=(taskset -c 0-3 nice -n 5)
# Waits (up to 30 min) until no other protal, simulator or test run is going (other sessions share the machine); prints the
# 1-minute load then.
quiet() { local n=0; while pgrep -f "build/protal |build/simulate_metagenomes|python3 -m unittest|build_gtdb_database|collect_training_data" > /dev/null && [ $n -lt 360 ]; do sleep 5; n=$((n + 1)); done; cut -d" " -f1 /proc/loadavg; }
secs() { awk '{ t = 0; for (i = 1; i <= NF; i++) { if ($i ~ /^[0-9]+m$/) t += 60 * $i; else if ($i ~ /^[0-9]+s$/) t += $i; else if ($i ~ /^[0-9]+ms$/) t += $i / 1000 } printf "%.3f", t }'; }
mem() { grep -m1 "^Memory after $1:" "$2" | sed -E 's/.*: ([0-9.]+) GB resident, peak ([0-9.]+) GB/\1\/\2/'; }
printf 'sample\tvariant\trep\twall_s\tmax_rss_mb\ttables_s\tlater_s\tlater_wait_s\tpreload_s\tindex_s\taligning_s\tprofiling_s\tafter_index\tafter_aligning\tafter_profiling\tload\n' > "$R/runs.tsv"
for rep in $(seq 1 "$reps"); do
  for sample in pe pb; do
    for variant in ref_text work_text work_bin; do
      case $variant in ref_text) bin=$HOME/impl/ref/build/protal db=$W/db_ref ;; work_text) bin=$HOME/impl/work/build/protal db=$W/db_ref ;;
                       work_bin) bin=$HOME/impl/work/build/protal db=$W/db_bin ;; esac
      case $sample in pe) reads=(-1 "$S/pe1m/reads/pe1m_1_R1.fq.gz" -2 "$S/pe1m/reads/pe1m_1_R2.fq.gz") ;; pb) reads=(-1 "$S/long/pb.fq.gz" --read_type pb) ;; esac
      out=$R/work; rm -rf "$out"
      load=$(quiet)
      log=$R/logs/$sample.$variant.$rep.log
      /usr/bin/time -f "%e %M" -o "$R/time.txt" "${pin[@]}" "$bin" --db "$db/database.protal" "${reads[@]}" -o "$out" -t 4 --no_qcmsa > "$log" 2>&1 \
          || echo "FAILED $sample $variant $rep"
      read -r wall rss < <(tail -1 "$R/time.txt")
      later=$(grep -m1 "^Loading the profiling's tables beside the alignment took" "$log")
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$variant" "$rep" "$wall" "$((rss / 1024))" \
        "$(grep -m1 '^Loading the taxonomy, models and tables took' "$log" | secs)" \
        "$(echo "$later" | grep -oE 'took [0-9.]+' | grep -oE '[0-9.]+$' || echo NA)" \
        "$(echo "$later" | grep -oE 'waited [0-9.]+' | grep -oE '[0-9.]+$' || echo NA)" \
        "$(grep -m1 '^Preload genomes took' "$log" | secs)" "$(grep -m1 '^Load Index took' "$log" | secs)" \
        "$(grep -m1 '^Aligning reads took' "$log" | secs)" "$(grep -m1 '^Profiling took' "$log" | secs)" \
        "$(mem 'loading the index' "$log")" "$(mem aligning "$log")" "$(mem profiling "$log")" "$load" >> "$R/runs.tsv"
      rm -rf "$out"
    done
  done
done
# The cohort: 24 samples of 50,000 pairs, strains on, no qcMSA.
printf 'variant\twall_s\tmax_rss_mb\tafter_aligning\tafter_profiling\tload\n' > "$R/cohort.tsv"
for variant in ref_text work_bin; do
  case $variant in ref_text) bin=$HOME/impl/ref/build/protal db=$W/db_ref ;; work_bin) bin=$HOME/impl/work/build/protal db=$W/db_bin ;; esac
  out=$R/cohort_$variant; rm -rf "$out"; mkdir -p "$out"
  load=$(quiet)
  { printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\tREAD_TYPE\n' "$out"
    for i in $(seq 1 24); do printf 'c%s\t%s\t%s\tc%s\tpe\n' "$i" "$S/cohort/reads/c_${i}_R1.fq.gz" "$S/cohort/reads/c_${i}_R2.fq.gz" "$i"; done; } > "$R/cohort_$variant.map"
  log=$R/logs/cohort.$variant.log
  /usr/bin/time -f "%e %M" -o "$R/time.txt" "${pin[@]}" "$bin" --db "$db/database.protal" --map "$R/cohort_$variant.map" -t 4 --no_qcmsa > "$log" 2>&1 \
      || echo "FAILED cohort $variant"
  read -r wall rss < <(tail -1 "$R/time.txt")
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$variant" "$wall" "$((rss / 1024))" "$(mem aligning "$log")" "$(mem profiling "$log")" "$load" >> "$R/cohort.tsv"
done
cat "$R/runs.tsv" "$R/cohort.tsv"
grep -h -E "^(Column weights|Strain alleles|Ancestry sites|Congener gaps):" "$R/logs/pe.work_bin.1.log" "$R/logs/pe.ref_text.1.log" | cut -c1-260
echo RUNS DONE
