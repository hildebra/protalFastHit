#!/usr/bin/env bash
# measure_performance.sh - time protal runs on a database and some read files: whole runs, protal's
# stage timers and, where perf works, CPU counters. For comparing machines, builds and databases.
#
#   bash scripts/measure_performance.sh OUT_DIR DB TYPE:R1[:R2] [TYPE:R1[:R2] ...]
#
# TYPE is protal's --read_type (pe, se, pb, ont); pe takes R1 and R2. Each sample runs REPEATS times
# on THREADS threads, one sample at a time; protal's outputs are removed after each run (KEEP=1 keeps
# them). Writes OUT_DIR/environment.txt (machine, protal, database), runs.tsv (one line per run),
# stages.tsv (each run's misc/<prefix>_runtime.tsv), logs/, and prints medians per sample.
# The first run reads the database from disk unless the page cache holds it, so its load_index_s
# is the cold load; later runs load from memory.
#
# Environment: THREADS (default: SLURM_CPUS_PER_TASK, else nproc), REPEATS (default 3),
#              PROTAL (default: build/protal), PROTAL_ARGS (more options, e.g. "--no_strains"),
#              PERF=0 (no perf stat), KEEP=1.
# On a cluster, ask for a whole node so other jobs do not share its memory bandwidth, e.g.
#   sbatch --exclusive -c 32 --mem 100G --wrap "bash scripts/measure_performance.sh perf_r226 \
#       r226_db pe:sample_R1.fq.gz:sample_R2.fq.gz ont:sample_ont.fq.gz"
set -uo pipefail

[ $# -ge 3 ] || { sed -n '2,19p' "$0" >&2; exit 1; }
out=$1; db=$2; shift 2
threads=${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}
repeats=${REPEATS:-3}
protal=${PROTAL:-build/protal}
[ -x "$protal" ] || { echo "protal binary not found at $protal (set PROTAL)" >&2; exit 1; }
[ -e "$db" ] || { echo "database not found: $db" >&2; exit 1; }
for spec in "$@"; do
    IFS=: read -r type r1 r2 <<< "$spec"
    case $type in pe) [ -n "$r2" ] || { echo "$spec: pe needs R1 and R2" >&2; exit 1; } ;;
                  se|pb|ont) ;; *) echo "$spec: TYPE is pe, se, pb or ont" >&2; exit 1 ;; esac
    for f in "$r1" ${r2:+"$r2"}; do [ -r "$f" ] || { echo "cannot read $f" >&2; exit 1; }; done
done
mkdir -p "$out/logs"

use_time=0; [ -x /usr/bin/time ] && use_time=1
use_perf=0
if [ "${PERF:-1}" = 1 ] && command -v perf > /dev/null && perf stat -x, -e instructions -o /dev/null true 2> /dev/null; then
    use_perf=1
fi

{
    echo "date: $(date -Is)"; echo "host: $(hostname)"; echo "job: ${SLURM_JOB_ID:-none}"
    echo "threads: $threads (CPUs this process may use: $(nproc))"
    echo "protal: $("$protal" --version 2>&1 | head -3 | tr '\n' ' ')"
    echo "database: $db, $(du -shL "$db" 2> /dev/null | cut -f1)"
    df -hT "$db" 2> /dev/null | tail -1 | awk '{ print "database file system: " $2 " (" $1 ")" }'
    echo "transparent huge pages: $(cat /sys/kernel/mm/transparent_hugepage/enabled 2> /dev/null)"
    echo "perf stat: $([ $use_perf = 1 ] && echo yes || echo "no (not installed or not allowed: perf_event_paranoid $(cat /proc/sys/kernel/perf_event_paranoid 2> /dev/null))")"
    echo "samples: $*"
    echo; lscpu 2> /dev/null | grep -E '^(Model name|Socket|Core|Thread|NUMA node\(s\)|L3|CPU max MHz)'
    echo; free -g 2> /dev/null
} > "$out/environment.txt"
cat "$out/environment.txt"; echo

# A stage's seconds from its "<stage> took 1m 2s 345ms" line of the log (NA if it did not run).
took() {
    grep -m1 "$1 took" "$2" | sed -e 's/.* took //' -e 's/ (.*//' -e 's/ mean over.*//' | awk '
        { s = 0; for (i = 1; i <= NF; i++) { v = $i; n = v + 0
              if (v ~ /^[0-9]+ms$/) s += n / 1000; else if (v ~ /^[0-9]+h$/) s += n * 3600
              else if (v ~ /^[0-9]+m$/) s += n * 60; else if (v ~ /^[0-9]+s$/) s += n }
          printf "%.3f", s; found = 1 }
        END { if (!found) printf "NA" }'
}
# The sum of a perf counter over the CPU types that report it (cpu_core and cpu_atom on hybrid CPUs).
counter() { awk -F, -v e="$1" '$3 ~ e && $1 ~ /^[0-9.]+$/ { s += $1; n++ } END { if (n) printf "%.0f", s; else printf "NA" }' "$2"; }

printf 'sample\ttype\trep\tthreads\twall_s\tuser_s\tsys_s\tmax_rss_gb\tmajor_faults\tload_index_s\taligning_s\tprofiling_s\tstrains_s\tinstructions\tcycles\tcache_misses\tipc\n' > "$out/runs.tsv"
printf 'sample\trep\tstage\tseconds\tthreads\tseconds_per_thread\n' > "$out/stages.tsv"
n=0
for spec in "$@"; do
    n=$((n + 1)); IFS=: read -r type r1 r2 <<< "$spec"; name=$type$n
    reads=(-1 "$r1"); [ -n "$r2" ] && reads+=(-2 "$r2")
    for rep in $(seq 1 "$repeats"); do
        run=$out/run.$name.$rep; log=$out/logs/$name.$rep.log; rm -rf "$run"
        cmd=("$protal" --db "$db" "${reads[@]}" --read_type "$type" -t "$threads" --prefix "$name" -o "$run" --no_qcmsa)
        # shellcheck disable=SC2206
        cmd+=(${PROTAL_ARGS:-})
        [ $use_perf = 1 ] && cmd=(perf stat -x, -o "$out/perf.tmp" -e cycles,instructions,cache-misses -- "${cmd[@]}")
        start=$(date +%s.%N)
        if [ $use_time = 1 ]; then /usr/bin/time -f '%e %U %S %M %F' -o "$out/time.tmp" "${cmd[@]}" > "$log" 2>&1
        else "${cmd[@]}" > "$log" 2>&1; fi
        status=$?
        wall=$(awk -v a="$start" -v b="$(date +%s.%N)" 'BEGIN { printf "%.2f", b - a }')
        [ $status = 0 ] || { echo "protal failed on $name (run $rep), see $log" >&2; tail -5 "$log" >&2; continue; }
        user=NA; sys=NA; rss=NA; faults=NA
        [ $use_time = 1 ] && read -r wall user sys rss faults < <(tail -1 "$out/time.tmp") && rss=$(awk -v k="$rss" 'BEGIN { printf "%.1f", k / 1048576 }')
        ins=NA; cyc=NA; miss=NA; ipc=NA
        if [ $use_perf = 1 ]; then
            ins=$(counter instructions "$out/perf.tmp"); cyc=$(counter cycles "$out/perf.tmp"); miss=$(counter cache-misses "$out/perf.tmp")
            [ "$ins" != NA ] && [ "$cyc" != NA ] && ipc=$(awk -v i="$ins" -v c="$cyc" 'BEGIN { printf "%.2f", i / c }')
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$type" "$rep" "$threads" \
            "$wall" "$user" "$sys" "$rss" "$faults" "$(took 'Load Index' "$log")" "$(took 'Aligning reads' "$log")" \
            "$(took 'Profiling' "$log")" "$(took 'Strain-level MSAs' "$log")" "$ins" "$cyc" "$miss" "$ipc" >> "$out/runs.tsv"
        [ -e "$run/misc/${name}_runtime.tsv" ] && awk -v s="$name" -v r="$rep" 'NR > 1 { print s "\t" r "\t" $0 }' \
            "$run/misc/${name}_runtime.tsv" >> "$out/stages.tsv"
        echo "$name run $rep: ${wall} s"
        [ "${KEEP:-0}" = 1 ] || rm -rf "$run"
    done
done
rm -f "$out/time.tmp" "$out/perf.tmp"

# Medians per sample over its runs, and its largest stages (seconds per thread, median). Plain awk (no gawk).
median='function median(a, n,   i, j, x) {
    for (i = 2; i <= n; i++) { x = a[i]; for (j = i - 1; j >= 1 && a[j] > x; j--) a[j + 1] = a[j]; a[j + 1] = x }
    return n % 2 ? a[(n + 1) / 2] : (a[n / 2] + a[n / 2 + 1]) / 2 }'
echo
awk -F'\t' "$median"'
    NR == 1 { for (i = 5; i <= 13; i++) h[i] = $i; next }
    { k = $1; if (!(k in c)) order[++m] = k; c[k]++; for (i = 5; i <= 13; i++) v[k, i, c[k]] = $i }
    END { for (j = 1; j <= m; j++) { k = order[j]; printf "%s (median of %d):", k, c[k]
            for (i = 5; i <= 13; i++) { split("", a); n = 0
                for (r = 1; r <= c[k]; r++) if (v[k, i, r] != "NA") a[++n] = v[k, i, r] + 0
                printf " %s %s", h[i], n ? sprintf("%.2f", median(a, n)) : "NA" }
            print "" } }' "$out/runs.tsv"
echo
awk -F'\t' "$median"'
    NR > 1 { k = $1 "\t" $3; if (!(k in n)) keys[++m] = k; v[k, ++n[k]] = $6 + 0 }
    END { for (j = 1; j <= m; j++) { k = keys[j]; split("", a); for (i = 1; i <= n[k]; i++) a[i] = v[k, i]
            printf "%s\t%.2f\n", k, median(a, n[k]) } }' "$out/stages.tsv" |
    sort -t$'\t' -k1,1 -k3,3gr | awk -F'\t' '{ if (++c[$1] <= 6) printf "%s  %-40s %8s s per thread\n", $1, $2, $3 }'
echo; echo "Results in $out (runs.tsv, stages.tsv, environment.txt, logs/)."
