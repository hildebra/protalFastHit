#!/usr/bin/env bash
# measure_performance.sh - time protal runs on a database and some read files: whole runs, protal's
# stage timers and, where perf works, CPU counters. For comparing machines, builds and databases.
#
#   bash scripts/measure_performance.sh [--threads N] [--sylph SYLPH_DB] OUT_DIR DB TYPE:R1[:R2] [TYPE:R1[:R2] ...]
#
# TYPE is protal's --read_type (pe, se, pb, ont); pe takes R1 and R2. Each sample runs REPEATS times
# on N threads, one sample at a time (--no_qcmsa: one sample makes no strain MSAs); protal's
# outputs are removed after each run (KEEP=1 keeps them). Then, with two samples or more, the cohort:
# all samples in one run from a map (COHORT_REPEATS times), which reuses the samples' SAM files, so
# that it times the profiling of all samples, the strain MSAs and qcMSA (QCMSA=0: --no_qcmsa).
# Writes OUT_DIR/environment.txt (machine, protal, database), runs.tsv (one line per sample run:
# times, memory, perf counters, protal's counts of reads, anchors and candidate alignments, the SAM
# header's genes and finishing time, the seeding's k-mers, blocks and flex cells, the seeds that share
# their taxon and gene with another, the lookups dropped as too ubiquitous, the anchors), the candidates the
# indel bound refused (ONT reads; NA for a protal from before 2026-10-09) and, last, the start-up's tables (tables_s:
# what the alignment waits for; profiling_tables_s, profiling_tables_wait_s: the tables loaded beside the alignment
# since 2026-10-11, and how long the profiling waited for them; NA before)), cohort.tsv (one
# line per cohort run: profiling, strain MSAs, building the MSAs, qcMSA, species), stages.tsv (each
# run's misc/<prefix>_runtime.tsv: the alignment stage's timers per thread and the profiling steps'
# wall times), logs/, and prints medians.
# The first run reads the database from disk unless the page cache holds it, so its load_index_s
# is the cold load; later runs load from memory.
# --threads N (or -t N): the threads protal and sylph each run on (their -t); without it, THREADS below.
# --sylph SYLPH_DB (a .syldb or .syl2db, best the same GTDB release as DB): last, sylph profiles each sample from
# the same reads, REPEATS times on N threads (sylph profile SYLPH_DB -1 R1 -2 R2, or -r R1 for se, pb and ont),
# for comparison: sylph.tsv (one line per sylph run: times, memory, perf counters, the genomes in sylph's profile),
# the profiles in sylph/, the logs in logs/sylph.*, and medians next to protal's. Its first run is the cold load too.
#
# Environment: THREADS (without --threads; default: SLURM_CPUS_PER_TASK, else nproc), REPEATS (default 3),
#              PROTAL (default: build/protal), PROTAL_ARGS (more options, e.g. "--no_strains"),
#              PERF=0 (no perf stat), KEEP=1, COHORT=0 (no cohort run), COHORT_REPEATS (default 1),
#              QCMSA=0 (the cohort without qcMSA; qcmsa must be next to protal, on $PATH or in
#              PROTAL_QCMSA_SCRIPT, see protal --help),
#              SYLPH (default: sylph on $PATH, e.g. from its activated conda environment),
#              SYLPH_ARGS (more sylph profile options, e.g. "-u").
# On a cluster, ask for a whole node so other jobs do not share its memory bandwidth, e.g.
#   sbatch --exclusive -c 32 --mem 100G --wrap "bash scripts/measure_performance.sh perf_r226 \
#       r226_db pe:sample_R1.fq.gz:sample_R2.fq.gz ont:sample_ont.fq.gz"
set -uo pipefail

# --threads N (-t N, --threads=N) and --sylph SYLPH_DB (--sylph=SYLPH_DB) anywhere among the arguments; the rest are
# positional. Without enough of them, the comment above.
threads=; sylph_db=; positional=()
while [ $# -gt 0 ]; do
    case $1 in
        -t|--threads) [ $# -ge 2 ] || { echo "$1 needs a number of threads" >&2; exit 1; }
                      threads=$2; shift 2 ;;
        --threads=*) threads=${1#--threads=}; shift ;;
        --sylph) [ $# -ge 2 ] || { echo "--sylph needs sylph's database (.syldb or .syl2db)" >&2; exit 1; }
                 sylph_db=$2; shift 2 ;;
        --sylph=*) sylph_db=${1#--sylph=}; shift ;;
        *) positional+=("$1"); shift ;;
    esac
done
set -- ${positional[@]+"${positional[@]}"}
[ $# -ge 3 ] || { awk 'NR > 1 { if (!/^#/) exit; print }' "$0" >&2; exit 1; }
out=$1; db=$2; shift 2
threads=${threads:-${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}}
[[ $threads =~ ^[1-9][0-9]*$ ]] || { echo "threads: a positive number, not '$threads'" >&2; exit 1; }
repeats=${REPEATS:-3}
cohort_repeats=${COHORT_REPEATS:-1}
protal=${PROTAL:-build/protal}
[ -x "$protal" ] || { echo "protal binary not found at $protal (set PROTAL)" >&2; exit 1; }
[ -e "$db" ] || { echo "database not found: $db" >&2; exit 1; }
for spec in "$@"; do
    IFS=: read -r type r1 r2 <<< "$spec"
    case $type in pe) [ -n "$r2" ] || { echo "$spec: pe needs R1 and R2" >&2; exit 1; } ;;
                  se|pb|ont) ;; *) echo "$spec: TYPE is pe, se, pb or ont" >&2; exit 1 ;; esac
    for f in "$r1" ${r2:+"$r2"}; do [ -r "$f" ] || { echo "cannot read $f" >&2; exit 1; }; done
done
sylph=; sylph_reads_option=0
if [ -n "$sylph_db" ]; then
    sylph=$(command -v "${SYLPH:-sylph}" || true)
    [ -n "$sylph" ] || { echo "sylph not found (activate its conda environment or set SYLPH)" >&2; exit 1; }
    [ -e "$sylph_db" ] || { echo "sylph database not found: $sylph_db" >&2; exit 1; }
    # Single-end reads with -r where sylph profile has it; an older sylph takes them as positional files (reads if FASTQ).
    case $("$sylph" profile --help 2>&1) in *--reads*) sylph_reads_option=1 ;; esac
fi
mkdir -p "$out/logs"
abs() { case $1 in /*) printf '%s' "$1" ;; *) printf '%s/%s' "$PWD" "$1" ;; esac; }

use_time=0; [ -x /usr/bin/time ] && use_time=1
use_perf=0
if [ "${PERF:-1}" = 1 ] && command -v perf > /dev/null && perf stat -x, -e instructions -o /dev/null true 2> /dev/null; then
    use_perf=1
fi
cohort=0; [ "${COHORT:-1}" = 1 ] && [ $# -ge 2 ] && cohort=1
qcmsa=${QCMSA:-1}

{
    echo "date: $(date -Is)"; echo "host: $(hostname)"; echo "job: ${SLURM_JOB_ID:-none}"
    echo "threads: $threads (CPUs this process may use: $(nproc))"
    echo "protal: $("$protal" --version 2>&1 | head -3 | tr '\n' ' ')"
    echo "database: $db, $(du -shL "$db" 2> /dev/null | cut -f1)"
    df -hT "$db" 2> /dev/null | tail -1 | awk '{ print "database file system: " $2 " (" $1 ")" }'
    echo "transparent huge pages: $(cat /sys/kernel/mm/transparent_hugepage/enabled 2> /dev/null)"
    echo "perf stat: $([ $use_perf = 1 ] && echo yes || echo "no (not installed or not allowed: perf_event_paranoid $(cat /proc/sys/kernel/perf_event_paranoid 2> /dev/null))")"
    echo "samples: $*"
    if [ $cohort = 1 ]; then
        # qcmsa as protal finds it (RunProtal.h, FindQCMSAScript): PROTAL_QCMSA_SCRIPT, qcmsa next to protal, qcmsa on
        # $PATH, then qcmsa.py next to protal, in its scripts/ or in the scripts/ beside its folder (a source checkout).
        bin=$(dirname "$protal"); q=${PROTAL_QCMSA_SCRIPT:-}
        [ -n "$q" ] || { [ -x "$bin/qcmsa" ] && q=$bin/qcmsa; }
        [ -n "$q" ] || q=$(command -v qcmsa || true)
        for c in "$bin/qcmsa.py" "$bin/scripts/qcmsa.py" "$bin/../scripts/qcmsa.py"; do [ -n "$q" ] || { [ -f "$c" ] && q=$c; }; done
        echo "cohort: all $# samples, $cohort_repeats run(s), qcMSA $([ "$qcmsa" = 1 ] && echo "on (${q:-qcmsa not found: protal will report it})" || echo off)"
    else
        echo "cohort: none ($([ $# -lt 2 ] && echo "one sample" || echo "COHORT=0"))"
    fi
    if [ -n "$sylph" ]; then
        echo "sylph: $("$sylph" --version 2>&1 | head -1) ($sylph), database $sylph_db, $(du -shL "$sylph_db" 2> /dev/null | cut -f1)"
    else
        echo "sylph: none"
    fi
    echo; lscpu 2> /dev/null | grep -E '^(Model name|Socket|Core|Thread|NUMA node\(s\)|L3|CPU max MHz)'
    echo; free -g 2> /dev/null
} > "$out/environment.txt"
cat "$out/environment.txt"; echo

# A stage's seconds from its "<stage> took 1m 2s 345ms" line of the log (NA if it did not run); anywhere on a line, as a
# message can follow qcmsa's progress bars, which end without a newline.
took() {
    grep -m1 -oE "$1 took .*" "$2" | sed -e "s/^$1 took //" -e 's/ (.*//' -e 's/ mean over.*//' | awk '
        { s = 0; for (i = 1; i <= NF; i++) { v = $i; n = v + 0
              if (v ~ /^[0-9]+ms$/) s += n / 1000; else if (v ~ /^[0-9]+h$/) s += n * 3600
              else if (v ~ /^[0-9]+m$/) s += n * 60; else if (v ~ /^[0-9]+s$/) s += n }
          printf "%.3f", s; found = 1 }
        END { if (!found) printf "NA" }'
}
# The sum of a perf counter over the CPU types that report it (cpu_core and cpu_atom on hybrid CPUs).
counter() { awk -F, -v e="$1" '$3 ~ e && $1 ~ /^[0-9.]+$/ { s += $1; n++ } END { if (n) printf "%.0f", s; else printf "NA" }' "$2"; }
# The first number after `pattern` (an extended regex) in the first line matching `line`, or NA.
number_after() {
    local v
    v=$(grep -m1 -E "$1" "$3" | grep -oE "$2 *[0-9.]+" | head -1 | grep -oE '[0-9.]+$')
    printf '%s' "${v:-NA}"
}
# The peak memory protal logged after a stage ("Memory after <stage>: R GB resident, peak P GB"), or NA.
peak_after() {
    number_after "^Memory after $1: " 'peak' "$2"
}
# The last timed run's perf counters (with use_perf=1, else NA): sets ins cyc miss ipc.
perf_counters() {
    ins=NA; cyc=NA; miss=NA; ipc=NA
    [ $use_perf = 1 ] || return 0
    ins=$(counter instructions "$out/perf.tmp"); cyc=$(counter cycles "$out/perf.tmp"); miss=$(counter cache-misses "$out/perf.tmp")
    [ "$ins" != NA ] && [ "$cyc" != NA ] && ipc=$(awk -v i="$ins" -v c="$cyc" 'BEGIN { printf "%.2f", i / c }')
    return 0
}
# Runs a command with its log, under /usr/bin/time (and perf stat with use_perf=1): sets wall user sys rss faults.
timed() {
    local log=$1 with_perf=$2; shift 2
    local cmd=("$@")
    [ "$with_perf" = 1 ] && cmd=(perf stat -x, -o "$out/perf.tmp" -e cycles,instructions,cache-misses -- "${cmd[@]}")
    local start; start=$(date +%s.%N)
    if [ $use_time = 1 ]; then /usr/bin/time -f '%e %U %S %M %F' -o "$out/time.tmp" "${cmd[@]}" > "$log" 2>&1
    else "${cmd[@]}" > "$log" 2>&1; fi
    status=$?
    wall=$(awk -v a="$start" -v b="$(date +%s.%N)" 'BEGIN { printf "%.2f", b - a }')
    user=NA; sys=NA; rss=NA; faults=NA
    [ $status = 0 ] && [ $use_time = 1 ] && read -r wall user sys rss faults < <(tail -1 "$out/time.tmp") &&
        rss=$(awk -v k="$rss" 'BEGIN { printf "%.1f", k / 1048576 }')
    return $status
}

printf 'sample\ttype\trep\tthreads\twall_s\tuser_s\tsys_s\tmax_rss_gb\tmajor_faults\tload_index_s\taligning_s\tprofiling_s\tstrains_s\tinstructions\tcycles\tcache_misses\tipc\treads\tanchored_reads\ttried\tscreened\tfrom_anchors\twhole_windows\tmade\twritten\tsam_finish_s\theader_genes\trecords_copy_s\tkmers\tkmers_in_index\tblocks_scanned\tflex_cells\tseeds\tpeak_preload_gb\tpeak_index_gb\tpeak_aligning_gb\tpeak_profiling_gb\tpaired_seeds\tdropped_lookups\tanchors\tindel_refused\ttables_s\tprofiling_tables_s\tprofiling_tables_wait_s\n' > "$out/runs.tsv"
printf 'sample\trep\tstage\tseconds\tthreads\tseconds_per_thread\n' > "$out/stages.tsv"
sams=$out/sams; mkdir -p "$sams"
cohort_map_rows=()
n=0
for spec in "$@"; do
    n=$((n + 1)); IFS=: read -r type r1 r2 <<< "$spec"; name=$type$n
    reads=(-1 "$r1"); [ -n "$r2" ] && reads+=(-2 "$r2")
    for rep in $(seq 1 "$repeats"); do
        run=$out/run.$name.$rep; log=$out/logs/$name.$rep.log; rm -rf "$run"
        cmd=("$protal" --db "$db" "${reads[@]}" --read_type "$type" -t "$threads" --prefix "$name" -o "$run" --no_qcmsa)
        # shellcheck disable=SC2206
        cmd+=(${PROTAL_ARGS:-})
        timed "$log" $use_perf "${cmd[@]}" || { echo "protal failed on $name (run $rep), see $log" >&2; tail -5 "$log" >&2; continue; }
        perf_counters
        # protal's counts line: "Sample <name>: R reads, A with an anchor; T candidate alignments tried: S refused by the k-mer
        # screen, I by the indel bound, F aligned from the anchor's exact matches and W as whole windows; M alignments made,
        # O records written" (a protal from before 2026-10-09 has no "I by the indel bound": its indel_refused is NA).
        numbers=$(grep -m1 "^Sample $name: .* candidate alignments tried" "$log" | sed 's/^Sample [^:]*: //' | grep -oE '[0-9]+' | paste -sd '\t' -)
        indel=NA
        case "$(printf '%s' "$numbers" | tr -cd '\t' | wc -c)" in
            8) indel=$(printf '%s' "$numbers" | cut -f5); counts=$(printf '%s' "$numbers" | cut -f1-4,6-9) ;;
            7) counts=$numbers ;;
            *) counts=$(printf 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA') ;;
        esac
        # The SAM header line ("SAM header: G genes, ...; the records (B bytes) were copied behind it in X s, ...").
        header_genes=$(number_after '^SAM header: ' 'SAM header:' "$log")
        copy_s=$(grep -m1 '^SAM header: ' "$log" | grep -oE 'copied behind it in [0-9.]+' | grep -oE '[0-9.]+$'); copy_s=${copy_s:-0}
        grep -q '^SAM header: ' "$log" || copy_s=NA
        # The seeding line ("Sample <name> seeding: K k-mers looked up, F in the index, S blocks scanned with C flex cells ...; N seeds,
        # P of them sharing their taxon and gene with another seed; D lookups dropped as too ubiquitous; A anchors"; a protal
        # from before 2026-10-06 ends it at "N seeds", and the last three columns are NA).
        seedline="^Sample $name seeding: "
        seeding=$(printf '%s\t%s\t%s\t%s\t%s' "$(number_after "$seedline" 'seeding:' "$log")" "$(number_after "$seedline" 'looked up,' "$log")" \
                  "$(number_after "$seedline" 'in the index,' "$log")" "$(number_after "$seedline" 'scanned with' "$log")" \
                  "$(grep -m1 -E "$seedline" "$log" | grep -oE '[0-9]+ seeds(,|$)' | grep -oE '^[0-9]+' || echo NA)")
        seed_fates=$(printf '%s\t%s\t%s' "$(number_after "$seedline" 'seeds,' "$log")" "$(number_after "$seedline" 'with another seed;' "$log")" \
                     "$(grep -m1 -E "$seedline" "$log" | grep -oE '[0-9]+ anchors$' | grep -oE '^[0-9]+' || echo NA)")
        # The start-up's tables: what the alignment waits for, and (since 2026-10-11) those loaded beside it ("Loading the
        # profiling's tables beside the alignment took X s; the profiling waited Y s for them"; NA without).
        later="^Loading the profiling's tables beside the alignment took "
        tables=$(printf '%s\t%s\t%s' "$(took 'Loading the taxonomy, models and tables' "$log")" "$(number_after "$later" 'took' "$log")" \
                 "$(number_after "$later" 'waited' "$log")")
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$type" "$rep" "$threads" \
            "$wall" "$user" "$sys" "$rss" "$faults" "$(took 'Load Index' "$log")" "$(took 'Aligning reads' "$log")" \
            "$(took 'Profiling' "$log")" "$(took 'Strain-level MSAs' "$log")" "$ins" "$cyc" "$miss" "$ipc" "$counts" \
            "$(took 'Writing the SAM header and file' "$log")" "$header_genes" "$copy_s" "$seeding" \
            "$(peak_after 'the genome preload' "$log")" "$(peak_after 'loading the index' "$log")" "$(peak_after aligning "$log")" \
            "$(peak_after profiling "$log")" "$seed_fates" "$indel" "$tables" >> "$out/runs.tsv"
        [ -e "$run/misc/${name}_runtime.tsv" ] && awk -v s="$name" -v r="$rep" 'NR > 1 { print s "\t" r "\t" $0 }' \
            "$run/misc/${name}_runtime.tsv" >> "$out/stages.tsv"
        echo "$name run $rep: ${wall} s"
        # The last run's SAM for the cohort run, which profiles it again with the other samples'.
        if [ $cohort = 1 ] && [ "$rep" = "$repeats" ]; then
            sam=$(ls "$run/$name".sam* 2> /dev/null | head -1)
            if [ -n "$sam" ]; then
                mv "$sam" "$sams/" && cohort_map_rows+=("$(printf '%s\t%s\t%s\t%s\t%s\t%s' "$name" "$(abs "$r1")" "$( [ -n "$r2" ] && abs "$r2" || printf -- '-')" \
                                                          "$(abs "$sams/$(basename "$sam")")" "$name" "$type")")
            fi
        fi
        [ "${KEEP:-0}" = 1 ] || rm -rf "$run"
    done
done

# The cohort: every sample in one run from a map, their SAMs given (protal profiles them without aligning), so that the
# strain MSAs (a species in two samples or more) and qcMSA run. Logs as cohort.<rep>.log, one line per run in cohort.tsv.
printf 'rep\tsamples\tthreads\twall_s\tuser_s\tsys_s\tmax_rss_gb\tload_index_s\tprofiling_s\tstrains_s\tbuild_msa_s\tqcmsa_s\tspecies\traw_msas\tfiltered_msas\n' > "$out/cohort.tsv"
if [ $cohort = 1 ] && [ ${#cohort_map_rows[@]} -lt $# ]; then
    echo "Cohort run skipped: only ${#cohort_map_rows[@]} of the $# samples have a SAM" >&2
    cohort=0
fi
if [ $cohort = 1 ]; then
    for rep in $(seq 1 "$cohort_repeats"); do
        run=$out/run.cohort.$rep; log=$out/logs/cohort.$rep.log; map=$out/cohort.$rep.map; rm -rf "$run"
        { printf '#OUTPUT_DIR\t%s\n' "$(abs "$run")"; printf '#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tREAD_TYPE\n'; printf '%s\n' "${cohort_map_rows[@]}"; } > "$map"
        cmd=("$protal" --db "$db" --map "$map" -t "$threads")
        [ "$qcmsa" = 1 ] || cmd+=(--no_qcmsa)
        # shellcheck disable=SC2206
        cmd+=(${PROTAL_ARGS:-})
        timed "$log" 0 "${cmd[@]}" || { echo "protal failed on the cohort (run $rep), see $log" >&2; tail -5 "$log" >&2; continue; }
        grep -q 'All alignments are present' "$log" || echo "Warning: the cohort run aligned reads again (see $log)" >&2
        msa_line='^Strain-level MSAs of '
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$rep" "$#" "$threads" "$wall" "$user" "$sys" "$rss" \
            "$(took 'Load Index' "$log")" "$(took 'Profiling' "$log")" "$(took 'Strain-level MSAs' "$log")" \
            "$(took 'Building the strain MSAs' "$log")" "$(took 'qcMSA' "$log")" "$(number_after "$msa_line" 'MSAs of' "$log")" \
            "$(number_after "$msa_line" 'species:' "$log")" "$(number_after "$msa_line" 'raw MSAs,' "$log")" >> "$out/cohort.tsv"
        for f in "$run"/misc/*_runtime.tsv; do
            [ -e "$f" ] || continue
            s=$(basename "$f" _runtime.tsv)
            awk -v s="cohort:$s" -v r="$rep" 'NR > 1 { print s "\t" r "\t" $0 }' "$f" >> "$out/stages.tsv"
        done
        echo "cohort run $rep: ${wall} s"
        [ "${KEEP:-0}" = 1 ] || rm -rf "$run"
    done
fi
[ "${KEEP:-0}" = 1 ] || rm -rf "$sams"

# sylph (--sylph) on the same reads, for comparison: each sample REPEATS times, one line per run in sylph.tsv, the profiles
# (small) kept in sylph/<sample>.<rep>.tsv.
if [ -n "$sylph" ]; then
    mkdir -p "$out/sylph"
    printf 'sample\ttype\trep\tthreads\twall_s\tuser_s\tsys_s\tmax_rss_gb\tmajor_faults\tinstructions\tcycles\tcache_misses\tipc\tgenomes\n' > "$out/sylph.tsv"
    n=0
    for spec in "$@"; do
        n=$((n + 1)); IFS=: read -r type r1 r2 <<< "$spec"; name=$type$n
        if [ -n "$r2" ]; then reads=(-1 "$r1" -2 "$r2"); elif [ $sylph_reads_option = 1 ]; then reads=(-r "$r1"); else reads=("$r1"); fi
        for rep in $(seq 1 "$repeats"); do
            log=$out/logs/sylph.$name.$rep.log; profile=$out/sylph/$name.$rep.tsv; rm -f "$profile"
            cmd=("$sylph" profile "$sylph_db" "${reads[@]}" -t "$threads" -o "$profile")
            # shellcheck disable=SC2206
            cmd+=(${SYLPH_ARGS:-})
            timed "$log" $use_perf "${cmd[@]}" || { echo "sylph failed on $name (run $rep), see $log" >&2; tail -5 "$log" >&2; continue; }
            perf_counters
            # sylph's profile: a header, then one line per genome it found in the sample.
            genomes=NA; [ -f "$profile" ] && genomes=$(awk 'NR > 1 && NF { c++ } END { print c + 0 }' "$profile")
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$type" "$rep" "$threads" \
                "$wall" "$user" "$sys" "$rss" "$faults" "$ins" "$cyc" "$miss" "$ipc" "$genomes" >> "$out/sylph.tsv"
            echo "sylph $name run $rep: ${wall} s"
        done
    done
fi
rm -f "$out/time.tmp" "$out/perf.tmp"

# Medians per sample over its runs, and its largest stages (seconds per thread, median). Plain awk (no gawk).
median='function median(a, n,   i, j, x) {
    for (i = 2; i <= n; i++) { x = a[i]; for (j = i - 1; j >= 1 && a[j] > x; j--) a[j + 1] = a[j]; a[j + 1] = x }
    return n % 2 ? a[(n + 1) / 2] : (a[n / 2] + a[n / 2 + 1]) / 2 }'
echo
awk -F'\t' "$median"'
    NR == 1 { for (i = 5; i <= NF; i++) h[i] = $i; last = NF; next }
    { k = $1; if (!(k in c)) order[++m] = k; c[k]++; for (i = 5; i <= last; i++) v[k, i, c[k]] = $i }
    END { for (j = 1; j <= m; j++) { k = order[j]; printf "%s (median of %d):", k, c[k]
            for (i = 5; i <= last; i++) { if ((i > 13 && i < 26) || (i > 26 && i < 42)) continue; split("", a); n = 0
                for (r = 1; r <= c[k]; r++) if (v[k, i, r] != "NA") a[++n] = v[k, i, r] + 0
                printf " %s %s", h[i], n ? sprintf("%.2f", median(a, n)) : "NA" }
            print "" } }' "$out/runs.tsv"
echo
# protal's counts (the same in every run of a sample): from the first run.
awk -F'\t' 'NR == 1 { for (i = 18; i <= NF; i++) { h[i] = $i; if ($i == "indel_refused") ir = i } next }
    !($1 in seen) { seen[$1] = 1; printf "%s counts:", $1; for (i = 18; i <= 25; i++) printf " %s %s", h[i], $i
                    if (ir) printf " %s %s", h[ir], $ir
                    printf "\n%s SAM header: genes %s, records copied in %s s; seeding:", $1, $27, $28
                    for (i = 29; i <= NF && i <= 41; i++) if (i != ir) printf " %s %s", h[i], $i
                    if ($31 > 0) printf " (%.1f flex cells per block)", $32 / $31
                    print "" }' "$out/runs.tsv"
if [ "$(wc -l < "$out/cohort.tsv")" -gt 1 ]; then
    echo
    awk -F'\t' "$median"'
        NR == 1 { nf = NF; for (i = 4; i <= NF; i++) h[i] = $i; next }
        { c++; for (i = 4; i <= NF; i++) v[i, c] = $i; s = $2 }
        END { printf "cohort of %d samples (median of %d):", s, c
              for (i = 4; i <= nf; i++) { split("", a); n = 0
                  for (r = 1; r <= c; r++) if (v[i, r] != "NA") a[++n] = v[i, r] + 0
                  printf " %s %s", h[i], n ? sprintf("%.2f", median(a, n)) : "NA" }
              print "" }' "$out/cohort.tsv"
fi
# sylph's medians per sample and protal's against them (CPU = user + sys).
if [ -f "$out/sylph.tsv" ] && [ "$(wc -l < "$out/sylph.tsv")" -gt 1 ]; then
    echo
    awk -F'\t' "$median"'
        function med(x, t, k,   a, i, j) { split("", a); j = 0
            for (i = 1; i <= runs[t, k]; i++) if (x[t, k, i] != "NA") a[++j] = x[t, k, i] + 0
            return j ? median(a, j) : "NA" }
        function num(v) { return v == "NA" ? v : sprintf("%.2f", v) }
        function ratio(p, s) { return p == "NA" || s == "NA" || s == 0 ? "NA" : sprintf("%.2fx", p / s) }
        FNR == 1 { t = FILENAME ~ /sylph\.tsv$/ ? "sylph" : "protal"; next }
        { k = $1; if (t == "sylph" && !(k in seen)) { seen[k] = 1; order[++m] = k }
          r = ++runs[t, k]; wall[t, k, r] = $5; rss[t, k, r] = $8
          cpu[t, k, r] = ($6 == "NA" || $7 == "NA") ? "NA" : $6 + $7
          if (t == "sylph") genomes[t, k, r] = $14 }
        END { for (j = 1; j <= m; j++) { k = order[j]
                  sw = med(wall, "sylph", k); sc = med(cpu, "sylph", k); sr = med(rss, "sylph", k)
                  printf "sylph %s (median of %d): wall_s %s cpu_s %s max_rss_gb %s genomes %s\n", k, runs["sylph", k],
                      num(sw), num(sc), num(sr), med(genomes, "sylph", k)
                  if (runs["protal", k]) printf "%s protal/sylph: wall %s, CPU %s, memory %s\n", k,
                      ratio(med(wall, "protal", k), sw), ratio(med(cpu, "protal", k), sc), ratio(med(rss, "protal", k), sr) } }' \
        "$out/runs.tsv" "$out/sylph.tsv"
fi
echo
awk -F'\t' "$median"'
    NR > 1 { k = $1 "\t" $3; if (!(k in n)) keys[++m] = k; v[k, ++n[k]] = $6 + 0 }
    END { for (j = 1; j <= m; j++) { k = keys[j]; split("", a); for (i = 1; i <= n[k]; i++) a[i] = v[k, i]
            printf "%s\t%.2f\n", k, median(a, n[k]) } }' "$out/stages.tsv" |
    sort -t$'\t' -k1,1 -k3,3gr | awk -F'\t' '{ if (++c[$1] <= 6) printf "%s  %-40s %8s s per thread\n", $1, $2, $3 }'
echo; echo "Results in $out (runs.tsv, cohort.tsv, stages.tsv,$([ -n "$sylph" ] && echo " sylph.tsv, sylph/,") environment.txt, logs/)."
