#!/usr/bin/env bash
# measure_performance.sh - time protal runs on a database and some read files: whole runs, protal's
# stage timers and, where perf works, CPU counters. For comparing machines, builds and databases.
#
#   bash scripts/measure_performance.sh [--threads N] [--cold_reads] [--sylph SYLPH_DB] [--stage1_probe] [--test_sparse] OUT_DIR DB TYPE:R1[:R2] [TYPE:R1[:R2] ...]
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
# since 2026-10-11, and how long the profiling waited for them; NA before) and the alignment's CPU seconds, the cores
# it kept busy and the seconds a thread spent getting reads (aligning_cpu_s, aligning_cores, reader_s_per_thread, from
# misc/cpu.tsv and misc/<prefix>_runtime.tsv: few cores and much reader time mean the input's inflating is the limit),
# and from protal's input line (since 2026-10-13; NA before) how each read file was read (input: e.g.
# gzip-members:22606:4t,... for gzip of 22606 members by members on 4 threads; bgzf, zstd-frames, gzip-pieces, gzip,
# zstd, plain; +members@B where a file read in pieces went on by members from byte B, +one-thread@B where members fell
# back to one thread) and the share of the alignment's time its threads spent getting reads (reader_share; a run at
# 0.25 or more is flagged on its progress line: its input held it back)), cohort.tsv (one
# line per cohort run: profiling, strain MSAs, building the MSAs, qcMSA, species), stages.tsv (each
# run's misc/<prefix>_runtime.tsv: the alignment stage's timers per thread and the profiling steps'
# wall times), logs/, and prints medians.
# The first run reads the database from disk unless the page cache holds it, so its load_index_s
# is the cold load; later runs load from memory.
# --threads N (or -t N): the threads protal and sylph each run on (their -t); without it, THREADS below.
# --cold_reads: the read files are dropped from the page cache before each protal run (dd iflag=nocache, which needs no
# root), so that every run reads its input from the disk or network file system, as a sample's first and usually only
# run does; the database stays cached. Without it, only a first run on the node may read cold.
# --sylph SYLPH_DB (a .syldb or .syl2db, best the same GTDB release as DB): last, sylph profiles each sample from
# the same reads, REPEATS times on N threads (sylph profile SYLPH_DB -1 R1 -2 R2, or -r R1 for se, pb and ont),
# for comparison: sylph.tsv (one line per sylph run: times, memory, perf counters, the genomes in sylph's profile),
# the profiles in sylph/, the logs in logs/sylph.*, and medians next to protal's. Its first run is the cold load too.
# --stage1_probe: after a sample's timed runs, one more run with protal's --stage1_probe (Stage1Probe.h: what the first
# stage of a two-stage protal would keep), not counted in runs.tsv: its tables in stage1_probe/<sample>_stage1_probe.tsv,
# _stage1_probe_stride.tsv (since 2026-10-13: the rules for a stride pick of the cores instead of a hash sample),
# _stage1_probe_taxa.tsv and (short reads) _first_lookups_probe.tsv, its log in logs/<sample>.probe.log, and a line in
# probe.tsv (times and memory, the probe's reads and reads with records, and whether its profile is byte-identical to the
# last timed run's, as the option promises; a run whose profile differs is kept for a look). Printed last: each sample's
# probe line and the probe run's wall time against the timed runs' median.
# --test_sparse: after a sample's timed runs, one run at each of four gates of protal's --sparse_first C (a read looked
# up by 1/C of its k-mer cores first and given up without a near-exact hit; ChainAnchorFinder::SetSparseFirst), by
# default with the grid pick on the gate's sparse index (since 2026-10-17; before, the database's pick on the whole
# index): levels 2g, 4g, 8g, 16g for short reads (SPARSE_SHORT), 4g, 8g, 16g, 32g for long reads (SPARSE_LONG); not
# counted in runs.tsv. A level C with a g gates with the grid pick (--sparse_first_pick grid: the cores of grid level
# log2 C or more, which the database chose per gene) and loads only those entries (--sparse_index auto), so its
# memory, load and lookup times are such an index's (docs/claude/2026-10-13-gene-grid-index). It needs an index with
# grid levels (built since 2026-10-13): where the first run's log says the index has none, the g levels are left out,
# with a warning (also in environment.txt and at the end). Other levels: one without a letter (1 to 32) gates with the
# database's pick on the whole index (stride where build_settings.tsv records none, as at r226 v24; SPARSE_SHORT="1 2 4
# 8" SPARSE_LONG="4 8 16 32" are the levels of before 2026-10-17); one with an s (e.g. 8s, since 2026-10-13) gives
# protal's --sparse_seeding with the hash pick: the reads kept are seeded from the gate's cores alone, as an index of
# only those would seed them, and the log tells what such an index would hold (docs/claude/2026-10-13-sparse-seeding);
# 0 is an ungated run. On a database that records gates (build_settings.tsv, since 2026-10-12) the timed runs follow
# them: add 0 to the levels for an ungated run. A line each in
# sparse.tsv: times, memory, perf counters, the reads given up and the sampled lookups, the records written against the
# last timed run's, its profile against that run's (species called by both, lost, gained, Bray-Curtis of the
# abundances as written and over the species alone, the two unknown shares), the alignment's CPU, cores and reader time as in runs.tsv, and the gate's own seconds a thread (protal's
# "Sparse first" timer, from 2026-10-12; NA before); each level's stage timers in sparse_stages.tsv (as stages.tsv, with
# the level for the rep); the profiles in sparse/<sample>.c<C>.profile beside sparse/<sample>.plain.profile, the logs in
# logs/<sample>.sparse<C>.log. Printed last: each level's wall, aligning and CPU time against the timed runs' medians.
#
# Environment: THREADS (without --threads; default: SLURM_CPUS_PER_TASK, else nproc), REPEATS (default 3),
#              PROTAL (default: build/protal), PROTAL_ARGS (more options, e.g. "--no_strains"),
#              PERF=0 (no perf stat), KEEP=1, COHORT=0 (no cohort run), COHORT_REPEATS (default 1),
#              QCMSA=0 (the cohort without qcMSA; qcmsa must be next to protal, on $PATH or in
#              PROTAL_QCMSA_SCRIPT, see protal --help),
#              SYLPH (default: sylph on $PATH, e.g. from its activated conda environment),
#              SYLPH_ARGS (more sylph profile options, e.g. "-u"),
#              SPARSE_SHORT, SPARSE_LONG (--test_sparse's levels, default "2g 4g 8g 16g" and "4g 8g 16g 32g"),
#              SPARSE_ARGS (more options of the --test_sparse runs, e.g. "--sparse_first_score 15").
# On a cluster, ask for a whole node so other jobs do not share its memory bandwidth, e.g.
#   sbatch --exclusive -c 32 --mem 100G --wrap "bash scripts/measure_performance.sh perf_r226 \
#       r226_db pe:sample_R1.fq.gz:sample_R2.fq.gz ont:sample_ont.fq.gz"
# RUN_DIR (environment): the folder the runs write into, protal's in RUN_DIR/protal (each run's -o folder, the samples'
# SAMs for the cohort) and sylph's in RUN_DIR/sylph (its profiles, copied to OUT_DIR/sylph at the end); by default
# OUT_DIR. Best a disk of the node itself (e.g. RUN_DIR=$TMPDIR/perf): writing and closing a SAM on a network file
# system took 1-7 s of a 32-thread run (docs/claude/2026-10-12-feature-cost-vs-f1.md). The tables and logs stay in OUT_DIR.
set -uo pipefail

# --threads N (-t N, --threads=N), --sylph SYLPH_DB (--sylph=SYLPH_DB), --stage1_probe and --test_sparse anywhere among
# the arguments; the rest are positional. Without enough of them, the comment above.
threads=; sylph_db=; stage1_probe=0; test_sparse=0; cold_reads=0; positional=()
while [ $# -gt 0 ]; do
    case $1 in
        --cold_reads|--cold-reads) cold_reads=1; shift ;;
        --stage1_probe|--stage1-probe) stage1_probe=1; shift ;;
        --test_sparse|--test-sparse) test_sparse=1; shift ;;
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
if [ $stage1_probe = 1 ]; then
    case $("$protal" --full_help 2>&1) in
        *--stage1_probe*) ;;
        *) echo "$protal has no --stage1_probe (a protal from before 2026-10-12)" >&2; exit 1 ;;
    esac
fi
sparse_short=${SPARSE_SHORT:-2g 4g 8g 16g}; sparse_long=${SPARSE_LONG:-4g 8g 16g 32g}
if [ $test_sparse = 1 ]; then
    case $("$protal" --full_help 2>&1) in
        *--sparse_first*) ;;
        *) echo "$protal has no --sparse_first (a protal from before 2026-10-12)" >&2; exit 1 ;;
    esac
    for c in $sparse_short $sparse_long; do
        case $c in 0|1|2|4|8|16|32|1s|2s|4s|8s|16s|32s|2g|4g|8g|16g|32g) ;; *) echo "--test_sparse levels are 0, 1, 2, 4, 8, 16 or 32, each with an optional s (--sparse_seeding; 2 to 32: or g, the grid pick on its sparse index), not '$c'" >&2; exit 1 ;; esac
        if [ "${c%g}" != "$c" ]; then
            case $("$protal" --full_help 2>&1) in
                *--sparse_index*) ;;
                *) echo "$protal has no --sparse_index (level $c; a protal from before 2026-10-13): give SPARSE_SHORT and SPARSE_LONG levels without a g, e.g. \"1 2 4 8\" and \"4 8 16 32\"" >&2; exit 1 ;;
            esac
        fi
        if [ "${c%s}" != "$c" ]; then
            case $("$protal" --full_help 2>&1) in
                *--sparse_seeding*) ;;
                *) echo "$protal has no --sparse_seeding (level $c; a protal from before 2026-10-13)" >&2; exit 1 ;;
            esac
        fi
    done
fi
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
runs=${RUN_DIR:-$out}; mkdir -p "$runs/protal" || { echo "cannot make $runs/protal (RUN_DIR)" >&2; exit 1; }
abs() { case $1 in /*) printf '%s' "$1" ;; *) printf '%s/%s' "$PWD" "$1" ;; esac; }
# A read file's format by its first bytes, as protal tells them apart (ThreadedGzStreambuf): BGZF (inflated on several
# threads), gzip, zstd (each on one thread) or plain.
input_format() {
    local h; h=$(head -c 16 "$1" 2> /dev/null | od -An -tx1 | tr -d ' \n')
    if [ "${h:0:8}" = 1f8b0804 ] && [ "${h:24:4}" = 4243 ]; then echo BGZF
    elif [ "${h:0:4}" = 1f8b ]; then echo gzip
    elif [ "${h:0:8}" = 28b52ffd ]; then echo zstd
    elif [ -n "$h" ]; then echo plain
    else echo unreadable; fi
}
# With --cold_reads: the read files' pages dropped from the page cache (posix_fadvise DONTNEED, by dd), so that the next
# run reads them from the disk.
drop_reads() {
    [ $cold_reads = 1 ] || return 0
    local f
    for f in "$@"; do dd if="$f" iflag=nocache count=0 status=none 2> /dev/null || echo "cannot drop $f from the page cache" >&2; done
}
# A run's input line ("Sample <name> input: <file> <how>; ...; the T alignment threads spent R s each getting reads, P% of
# the alignment's A s[: the input held the alignment back]", since 2026-10-13): each file's reading in short, files joined
# by "," (e.g. gzip-members:22606:4t; +members@B, +one-thread@B), a tab, and the reader share ("0.67"); NA NA without it.
input_of() {
    local line; line=$(grep -m1 "^Sample $1 input: " "$2")
    [ -n "$line" ] || { printf 'NA\tNA'; return 0; }
    printf '%s\n' "$line" | awk '{
        sub(/^Sample [^ ]* input: /, ""); n = split($0, part, "; "); out = ""; share = "NA"
        for (i = 1; i <= n; i++) {
            p = part[i]
            if (p ~ /^the [0-9]+ alignment threads spent/) {
                if (match(p, /[0-9.]+% of/)) share = sprintf("%.2f", substr(p, RSTART, RLENGTH - 4) / 100)
                continue
            }
            how = p; sub(/^[^ ]+ /, "", how)
            count = 0; if (match(how, /[0-9]+ (block|frame|member)/)) { count = substr(how, RSTART, RLENGTH); sub(/ .*/, "", count) }
            t = 1; if (match(how, / on [0-9]+ threads/)) t = substr(how, RSTART + 4, RLENGTH - 12)
            code = how ~ /^BGZF/ ? "bgzf" : how ~ /^zstd, .* decompressed on [0-9]+ threads/ ? "zstd-frames" : \
                   how ~ /by members on/ ? "gzip-members" : how ~ /in pieces/ ? "gzip-pieces" : how ~ /^gzip/ ? "gzip" : \
                   how ~ /^zstd/ ? "zstd" : "plain"
            code = code ":" count ":" t "t"
            if (match(how, /by members from byte [0-9]+/)) code = code "+members@" substr(how, RSTART + 21, RLENGTH - 21)
            if (match(how, /on one thread from byte [0-9]+/)) code = code "+one-thread@" substr(how, RSTART + 24, RLENGTH - 24)
            out = out (out == "" ? "" : ",") code
        }
        printf "%s\t%s", (out == "" ? "NA" : out), share }'
}

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
    echo "cold reads: $([ $cold_reads = 1 ] && echo "yes, the read files dropped from the page cache before each protal run" || echo no)"
    echo "protal: $("$protal" --version 2>&1 | head -3 | tr '\n' ' ')"
    echo "database: $db, $(du -shL "$db" 2> /dev/null | cut -f1)"
    df -hT "$db" 2> /dev/null | tail -1 | awk '{ print "database file system: " $2 " (" $1 ")" }'
    df -hT "$runs" 2> /dev/null | tail -1 | awk -v r="$runs" '{ print "run folder (RUN_DIR): " r ", protal/ and sylph/, file system " $2 " (" $1 ")" }'
    echo "transparent huge pages: $(cat /sys/kernel/mm/transparent_hugepage/enabled 2> /dev/null)"
    echo "perf stat: $([ $use_perf = 1 ] && echo yes || echo "no (not installed or not allowed: perf_event_paranoid $(cat /proc/sys/kernel/perf_event_paranoid 2> /dev/null))")"
    echo "samples: $*"
    printf 'read files:'; for spec in "$@"; do IFS=: read -r _ f1 f2 <<< "$spec"
        printf ' %s %s' "$(basename "$f1")" "$(input_format "$f1")"; [ -n "$f2" ] && printf ', %s %s' "$(basename "$f2")" "$(input_format "$f2")"; printf ';'
    done; echo
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
    echo "stage-1 probe: $([ $stage1_probe = 1 ] && echo "one more run per sample with --stage1_probe" || echo none)"
    echo "sparse first: $([ $test_sparse = 1 ] && echo "one run per level, --sparse_first $sparse_short (short reads), $sparse_long (long reads)${SPARSE_ARGS:+ with $SPARSE_ARGS} (g: the grid pick on its sparse index; s: --sparse_seeding, hash pick; no letter: the database's pick, whole index)" || echo none)"
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
# A run's alignment, from the misc tables of its output folder $1 and sample $2: "CPU seconds<TAB>cores busy<TAB>Sequence
# reader seconds per thread" (misc/cpu.tsv's "aligning <sample>" row, its CPU over its wall time; misc/<sample>_runtime.tsv;
# NA where missing). Few cores busy and much reader time: the threads wait for the input's inflating thread(s).
aligning_cpu() {
    local cpu_cores=$'NA\tNA' reader=NA
    [ -e "$1/misc/cpu.tsv" ] && cpu_cores=$(awk -F'\t' '$1 ~ /^aligning / { c = $3; w = $2 }
        END { if (c == "") { printf "NA\tNA"; exit } cores = w > 0 ? c / w : 0; printf "%.1f\t%.1f", c, cores }' "$1/misc/cpu.tsv")
    [ -e "$1/misc/$2_runtime.tsv" ] && reader=$(awk -F'\t' '$1 == "Sequence reader" { v = $4 }
        END { if (v == "") print "NA"; else printf "%.2f\n", v }' "$1/misc/$2_runtime.tsv")
    printf '%s\t%s' "$cpu_cores" "$reader"
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
# The last timed run's peak memory for its progress line, after its wall time: " 36.6Gb" (max_rss_gb, GiB); nothing
# without /usr/bin/time.
mem() { [ "$rss" = NA ] || printf ' %sGb' "$rss"; }

printf 'sample\ttype\trep\tthreads\twall_s\tuser_s\tsys_s\tmax_rss_gb\tmajor_faults\tload_index_s\taligning_s\tprofiling_s\tstrains_s\tinstructions\tcycles\tcache_misses\tipc\treads\tanchored_reads\ttried\tscreened\tfrom_anchors\twhole_windows\tmade\twritten\tsam_finish_s\theader_genes\trecords_copy_s\tkmers\tkmers_in_index\tblocks_scanned\tflex_cells\tseeds\tpeak_preload_gb\tpeak_index_gb\tpeak_aligning_gb\tpeak_profiling_gb\tpaired_seeds\tdropped_lookups\tanchors\tindel_refused\ttables_s\tprofiling_tables_s\tprofiling_tables_wait_s\taligning_cpu_s\taligning_cores\treader_s_per_thread\tinput\treader_share\n' > "$out/runs.tsv"
printf 'sample\trep\tstage\tseconds\tthreads\tseconds_per_thread\n' > "$out/stages.tsv"
# The tables of the steps this run leaves out, left in OUT_DIR by an earlier run there: removed, so that neither the
# summary below nor a reader takes them for this run's (a run without --sylph printed the sylph medians of an older one).
[ -n "$sylph" ] || rm -rf "$out/sylph.tsv" "$out/sylph" "$out"/logs/sylph.*.log
[ $stage1_probe = 1 ] || rm -rf "$out/probe.tsv" "$out/stage1_probe" "$out"/logs/*.probe.log
[ $test_sparse = 1 ] || rm -rf "$out/sparse.tsv" "$out/sparse_stages.tsv" "$out/sparse" "$out"/logs/*.sparse*.log
if [ $stage1_probe = 1 ]; then
    mkdir -p "$out/stage1_probe"
    printf 'sample\ttype\tthreads\twall_s\tuser_s\tsys_s\tmax_rss_gb\taligning_s\tprobe_reads\tprobe_reads_with_records\tprofile_identical\n' > "$out/probe.tsv"
fi
if [ $test_sparse = 1 ]; then
    mkdir -p "$out/sparse"
    printf 'sample\ttype\tthreads\tsparse_first\twall_s\tuser_s\tsys_s\tmax_rss_gb\taligning_s\tinstructions\tcycles\tipc\treads\tgiven_up\tgiven_up_share\tsampled_lookups\tkmers_kept_reads\twritten\twritten_plain\tspecies_plain\tspecies\tspecies_both\tspecies_lost\tspecies_gained\tbray_curtis\tbray_curtis_species\tunknown_plain\tunknown\taligning_cpu_s\taligning_cores\treader_s_per_thread\tgate_s_per_thread\n' > "$out/sparse.tsv"
    printf 'sample\tsparse_first\tstage\tseconds\tthreads\tseconds_per_thread\n' > "$out/sparse_stages.tsv"
fi
# Two profiles' species (their lineage column) and abundances, tab-separated: "plain sparse both lost gained bray_curtis
# bray_curtis_species unknown_plain unknown": the species of each, of both, lost and gained, Bray-Curtis of the abundances
# as written, the same over the species alone (each profile's species scaled to sum 1: a change of the unknown share
# alone moves the first, not this), and each profile's unknown share (its "?" line; NA without). Lines of 3 fields or
# more, comments left out. Membership is tested before a[] is read: reading an awk array element creates it (every
# species counted as in both before 2026-10-12's fix).
compare_profiles() {
    awk -F'\t' 'FNR == 1 { f++ } NF >= 3 && $1 == "?" { u[f] = $3 + 0; next }
        NF >= 3 && $1 !~ /^#/ { a[f, $2] = $3 + 0; seen[$2]; n[f]++; sum[f] += $3 }
        END { for (s in seen) { in1 = ((1, s) in a); in2 = ((2, s) in a); x = in1 ? a[1, s] : 0; y = in2 ? a[2, s] : 0
                  d += (x > y ? x - y : y - x); t += x + y
                  xs = sum[1] > 0 ? x / sum[1] : 0; ys = sum[2] > 0 ? y / sum[2] : 0; ds += (xs > ys ? xs - ys : ys - xs)
                  if (in1 && in2) both++; else if (in1) lost++; else gained++ }
              # (the comparisons before printf: a ">" among its arguments would redirect its output)
              bc = t > 0 ? sprintf("%.4f", d / t) : "NA"; bcs = (sum[1] > 0 && sum[2] > 0) ? sprintf("%.4f", ds / 2) : "NA"
              u1 = (1 in u) ? sprintf("%.4f", u[1]) : "NA"; u2 = (2 in u) ? sprintf("%.4f", u[2]) : "NA"
              printf "%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\n", n[1], n[2], both, lost, gained, bc, bcs, u1, u2 }' "$1" "$2"
}
# --test_sparse's g levels need an index with grid levels, which protal's log tells as it loads it ("Index features: ...,
# no grid levels (built before 2026-10-13)"); protal only reads them after its gene tables, so the script cannot ask
# before its first run. On the log of the first run that loaded the index: without grid levels, the g levels are left
# out (each run would stop, exit 2), with a warning here, in environment.txt and at the end (grid_dropped).
grid_checked=0; grid_dropped=
check_grid_levels() {
    { [ $test_sparse = 1 ] && [ $grid_checked = 0 ] && grep -q '^Index features: ' "$1"; } || return 0
    grid_checked=1
    grep -q '^Index features: .*no grid levels' "$1" || return 0
    local c short= long= dropped=0
    for c in $sparse_short; do if [ "${c%g}" != "$c" ]; then dropped=1; else short+=" $c"; fi; done
    for c in $sparse_long; do if [ "${c%g}" != "$c" ]; then dropped=1; else long+=" $c"; fi; done
    [ $dropped = 1 ] || return 0
    short=${short# }; long=${long# }
    grid_dropped="--test_sparse's levels with a g are left out: the index of $db has no grid levels (built before 2026-10-13). Levels without a letter gate with the database's pick on the whole index (SPARSE_SHORT=\"1 2 4 8\" SPARSE_LONG=\"4 8 16 32\"). Kept: short reads '${short:-none}', long reads '${long:-none}'."
    sparse_short=$short; sparse_long=$long
    echo "Warning: $grid_dropped" >&2
    echo "sparse first: $grid_dropped" >> "$out/environment.txt"
}
sams=$runs/protal/sams; mkdir -p "$sams"
cohort_map_rows=()
n=0
for spec in "$@"; do
    n=$((n + 1)); IFS=: read -r type r1 r2 <<< "$spec"; name=$type$n
    reads=(-1 "$r1"); [ -n "$r2" ] && reads+=(-2 "$r2")
    for rep in $(seq 1 "$repeats"); do
        run=$runs/protal/run.$name.$rep; log=$out/logs/$name.$rep.log; rm -rf "$run"
        cmd=("$protal" --db "$db" "${reads[@]}" --read_type "$type" -t "$threads" --prefix "$name" -o "$run" --no_qcmsa)
        # shellcheck disable=SC2206
        cmd+=(${PROTAL_ARGS:-})
        drop_reads "$r1" ${r2:+"$r2"}
        timed "$log" $use_perf "${cmd[@]}" || { echo "protal failed on $name (run $rep), see $log" >&2; tail -5 "$log" >&2; continue; }
        perf_counters
        check_grid_levels "$log"
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
        input=$(input_of "$name" "$log")
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$type" "$rep" "$threads" \
            "$wall" "$user" "$sys" "$rss" "$faults" "$(took 'Load Index' "$log")" "$(took 'Aligning reads' "$log")" \
            "$(took 'Profiling' "$log")" "$(took 'Strain-level MSAs' "$log")" "$ins" "$cyc" "$miss" "$ipc" "$counts" \
            "$(took 'Writing the SAM header and file' "$log")" "$header_genes" "$copy_s" "$seeding" \
            "$(peak_after 'the genome preload' "$log")" "$(peak_after 'loading the index' "$log")" "$(peak_after aligning "$log")" \
            "$(peak_after profiling "$log")" "$seed_fates" "$indel" "$tables" "$(aligning_cpu "$run" "$name")" "$input" >> "$out/runs.tsv"
        [ -e "$run/misc/${name}_runtime.tsv" ] && awk -v s="$name" -v r="$rep" 'NR > 1 { print s "\t" r "\t" $0 }' \
            "$run/misc/${name}_runtime.tsv" >> "$out/stages.tsv"
        # A run whose threads spent a quarter or more of the alignment getting reads: its input held it back.
        held=$(printf '%s' "$input" | awk -F'\t' '$2 != "NA" && $2 >= 0.25 { printf " (the input held the alignment back: %d%% of it getting reads; %s)", 100 * $2, $1 }')
        echo "$name run $rep: ${wall} s$(mem)$held"
        # The last run's SAM for the cohort run, which profiles it again with the other samples'.
        if [ $cohort = 1 ] && [ "$rep" = "$repeats" ]; then
            sam=$(ls "$run/$name".sam* 2> /dev/null | head -1)
            if [ -n "$sam" ]; then
                mv "$sam" "$sams/" && cohort_map_rows+=("$(printf '%s\t%s\t%s\t%s\t%s\t%s' "$name" "$(abs "$r1")" "$( [ -n "$r2" ] && abs "$r2" || printf -- '-')" \
                                                          "$(abs "$sams/$(basename "$sam")")" "$name" "$type")")
            fi
        fi
        # The last timed run's profile, which the --stage1_probe run's must equal and the --test_sparse runs' are compared with.
        [ $stage1_probe = 1 ] && [ "$rep" = "$repeats" ] && [ -e "$run/$name.profile" ] && cp "$run/$name.profile" "$out/stage1_probe/$name.plain.profile"
        [ $test_sparse = 1 ] && [ "$rep" = "$repeats" ] && [ -e "$run/$name.profile" ] && cp "$run/$name.profile" "$out/sparse/$name.plain.profile"
        [ "${KEEP:-0}" = 1 ] || rm -rf "$run"
    done
    # --stage1_probe: one more run, with the probe (outside runs.tsv); its tables, a line in probe.tsv, and its profile
    # compared with the timed runs'.
    if [ $stage1_probe = 1 ]; then
        run=$runs/protal/run.$name.probe; log=$out/logs/$name.probe.log; rm -rf "$run"
        cmd=("$protal" --db "$db" "${reads[@]}" --read_type "$type" -t "$threads" --prefix "$name" -o "$run" --no_qcmsa --stage1_probe)
        # shellcheck disable=SC2206
        cmd+=(${PROTAL_ARGS:-})
        if timed "$log" 0 "${cmd[@]}"; then
            for table in stage1_probe stage1_probe_stride stage1_probe_taxa first_lookups_probe; do
                [ -e "$run/misc/${name}_$table.tsv" ] && cp "$run/misc/${name}_$table.tsv" "$out/stage1_probe/"
            done
            [ -e "$out/stage1_probe/${name}_stage1_probe.tsv" ] || echo "Warning: the --stage1_probe run of $name wrote no misc/${name}_stage1_probe.tsv (see $log)" >&2
            identical=NA; plain=$out/stage1_probe/$name.plain.profile
            if [ -e "$plain" ] && [ -e "$run/$name.profile" ]; then
                if cmp -s "$plain" "$run/$name.profile"; then identical=yes; else identical=no; fi
            fi
            # "Sample <name> stage-1 probe: N reads (pairs as one), A with records; <rules>; ..."
            probeline="^Sample $name stage-1 probe: "
            with_records=$(grep -m1 -E "$probeline" "$log" | grep -oE '[0-9]+ with records' | grep -oE '^[0-9]+' || echo NA)
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$type" "$threads" "$wall" "$user" "$sys" "$rss" \
                "$(took 'Aligning reads' "$log")" "$(number_after "$probeline" 'probe:' "$log")" "${with_records:-NA}" "$identical" >> "$out/probe.tsv"
            echo "$name stage-1 probe run: ${wall} s$(mem), profile identical to the timed run's: $identical"
            if [ "$identical" = no ]; then
                echo "Warning: $name's profile with --stage1_probe differs from the timed run's: kept in $run" >&2
            elif [ "${KEEP:-0}" != 1 ]; then
                rm -rf "$run"
            fi
        else
            echo "protal --stage1_probe failed on $name, see $log" >&2; tail -5 "$log" >&2
            [ "${KEEP:-0}" = 1 ] || rm -rf "$run"
        fi
    fi
    # --test_sparse: one run at each level of --sparse_first (outside runs.tsv); a line each in sparse.tsv, the profile
    # compared with the last timed run's.
    if [ $test_sparse = 1 ]; then
        case $type in pb|ont) levels=$sparse_long ;; *) levels=$sparse_short ;; esac
        written_plain=$(awk -F'\t' -v s="$name" '$1 == s { w = $25 } END { print w == "" ? "NA" : w }' "$out/runs.tsv")
        for c in $levels; do
            run=$runs/protal/run.$name.sparse$c; log=$out/logs/$name.sparse$c.log; rm -rf "$run"
            cmd=("$protal" --db "$db" "${reads[@]}" --read_type "$type" -t "$threads" --prefix "$name" -o "$run" --no_qcmsa --sparse_first "${c%[sg]}")
            # A level C with an s: the reads kept seeded from the gate's cores alone (--sparse_seeding), as an index of
            # only those would seed them. Such an index holds the hash pick's cores: a protal that picks by stride by
            # default (since 2026-10-13) is told the hash pick.
            if [ "${c%s}" != "$c" ]; then
                cmd+=(--sparse_seeding)
                case $("$protal" --full_help 2>&1) in *--sparse_first_pick*) cmd+=(--sparse_first_pick hash) ;; esac
            fi
            # A level C with a g: the grid pick, on the sparse index of its gate (--sparse_index auto, protal's default).
            [ "${c%g}" != "$c" ] && cmd+=(--sparse_first_pick grid)
            # shellcheck disable=SC2206
            cmd+=(${SPARSE_ARGS:-} ${PROTAL_ARGS:-})
            if ! timed "$log" $use_perf "${cmd[@]}"; then
                echo "protal --sparse_first $c failed on $name, see $log" >&2; tail -5 "$log" >&2
                [ "${KEEP:-0}" = 1 ] || rm -rf "$run"
                continue
            fi
            perf_counters
            # "Sample <name> sparse first (1/C of the cores, score T, M hit(s)): G of R pairs (or reads) given up (p%), L sampled
            # lookups (l a pair); the seeding line counts the kept ones'"
            sparseline="^Sample $name sparse first "
            given=$(number_after "$sparseline" 'hit\(s\)\):' "$log")
            total=$(number_after "$sparseline" 'of' "$log")
            sampled=$(grep -m1 -E "$sparseline" "$log" | grep -oE '[0-9]+ sampled lookups' | grep -oE '^[0-9]+' || echo NA)
            share=$(awk -v g="$given" -v t="$total" 'BEGIN { print (g == "NA" || t == "NA" || t == 0) ? "NA" : sprintf("%.4f", g / t) }')
            written=$(grep -m1 "^Sample $name: .* candidate alignments tried" "$log" | grep -oE '[0-9]+ records written' | grep -oE '^[0-9]+' || echo NA)
            species=$(printf 'NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA')
            if [ -e "$run/$name.profile" ]; then
                cp "$run/$name.profile" "$out/sparse/$name.c$c.profile"
                [ -e "$out/sparse/$name.plain.profile" ] && species=$(compare_profiles "$out/sparse/$name.plain.profile" "$run/$name.profile")
            fi
            # The run's stage timers (misc/<sample>_runtime.tsv; "Sparse first" is the gate's own time, within "Seed- and
            # Anchor-finding" for short reads and "Alignment handler" for long reads).
            gate=NA
            if [ -e "$run/misc/${name}_runtime.tsv" ]; then
                awk -v s="$name" -v c="$c" 'NR > 1 { print s "\t" c "\t" $0 }' "$run/misc/${name}_runtime.tsv" >> "$out/sparse_stages.tsv"
                gate=$(awk -F'\t' '$1 == "Sparse first" { v = $4 } END { if (v == "") print "NA"; else printf "%.2f\n", v }' "$run/misc/${name}_runtime.tsv")
            fi
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$type" "$threads" "$c" \
                "$wall" "$user" "$sys" "$rss" "$(took 'Aligning reads' "$log")" "$ins" "$cyc" "$ipc" "${total:-NA}" "${given:-NA}" "$share" \
                "${sampled:-NA}" "$(number_after "^Sample $name seeding: " 'seeding:' "$log")" "${written:-NA}" "$written_plain" "$species" \
                "$(aligning_cpu "$run" "$name")" "$gate" >> "$out/sparse.tsv"
            echo "$name --sparse_first $c: ${wall} s$(mem), given up $given of $total, records $written (timed runs: $written_plain)"
            [ "${KEEP:-0}" = 1 ] || rm -rf "$run"
        done
    fi
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
        run=$runs/protal/run.cohort.$rep; log=$out/logs/cohort.$rep.log; map=$out/cohort.$rep.map; rm -rf "$run"
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
        echo "cohort run $rep: ${wall} s$(mem)"
        [ "${KEEP:-0}" = 1 ] || rm -rf "$run"
    done
fi
[ "${KEEP:-0}" = 1 ] || rm -rf "$sams"

# sylph (--sylph) on the same reads, for comparison: each sample REPEATS times, one line per run in sylph.tsv, the profiles
# (small) kept in sylph/<sample>.<rep>.tsv.
if [ -n "$sylph" ]; then
    mkdir -p "$runs/sylph"
    printf 'sample\ttype\trep\tthreads\twall_s\tuser_s\tsys_s\tmax_rss_gb\tmajor_faults\tinstructions\tcycles\tcache_misses\tipc\tgenomes\n' > "$out/sylph.tsv"
    n=0
    for spec in "$@"; do
        n=$((n + 1)); IFS=: read -r type r1 r2 <<< "$spec"; name=$type$n
        if [ -n "$r2" ]; then reads=(-1 "$r1" -2 "$r2"); elif [ $sylph_reads_option = 1 ]; then reads=(-r "$r1"); else reads=("$r1"); fi
        for rep in $(seq 1 "$repeats"); do
            log=$out/logs/sylph.$name.$rep.log; profile=$runs/sylph/$name.$rep.tsv; rm -f "$profile"
            cmd=("$sylph" profile "$sylph_db" "${reads[@]}" -t "$threads" -o "$profile")
            # shellcheck disable=SC2206
            cmd+=(${SYLPH_ARGS:-})
            timed "$log" $use_perf "${cmd[@]}" || { echo "sylph failed on $name (run $rep), see $log" >&2; tail -5 "$log" >&2; continue; }
            perf_counters
            # sylph's profile: a header, then one line per genome it found in the sample.
            genomes=NA; [ -f "$profile" ] && genomes=$(awk 'NR > 1 && NF { c++ } END { print c + 0 }' "$profile")
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$name" "$type" "$rep" "$threads" \
                "$wall" "$user" "$sys" "$rss" "$faults" "$ins" "$cyc" "$miss" "$ipc" "$genomes" >> "$out/sylph.tsv"
            echo "sylph $name run $rep: ${wall} s$(mem)"
        done
    done
    # With RUN_DIR elsewhere, the profiles (small) are kept with the tables.
    if [ "$(abs "$runs")" != "$(abs "$out")" ]; then mkdir -p "$out/sylph" && cp "$runs"/sylph/*.tsv "$out/sylph/" 2> /dev/null; fi
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
if [ -n "$sylph" ] && [ -f "$out/sylph.tsv" ] && [ "$(wc -l < "$out/sylph.tsv")" -gt 1 ]; then
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
# --stage1_probe: each sample's probe line, and the probe run against the timed runs' median wall time.
if [ $stage1_probe = 1 ] && [ "$(wc -l < "$out/probe.tsv")" -gt 1 ]; then
    echo
    awk -F'\t' "$median"'
        FNR == 1 { t = FILENAME ~ /probe\.tsv$/ ? "probe" : "runs"; next }
        t == "runs" { k = $1; wall[k, ++c[k]] = $5; next }
        { k = $1; split("", a); n = 0
          for (r = 1; r <= c[k]; r++) if (wall[k, r] != "NA") a[++n] = wall[k, r] + 0
          m = n ? median(a, n) : 0
          # (a ">" among printf'"'"'s arguments would redirect its output: the comparisons come first)
          timed = n ? sprintf("%.2f", m) : "NA"; more = ""
          if (m > 0 && $4 != "NA") more = sprintf(", %+.1f%%", 100 * ($4 / m - 1))
          printf "%s stage-1 probe run: wall %s s against %s s (median of the timed runs%s), max_rss_gb %s, profile identical: %s\n",
              k, $4, timed, more, $7, $11 }' \
        "$out/runs.tsv" "$out/probe.tsv"
    for log in "$out"/logs/*.probe.log; do
        [ -e "$log" ] && grep -m1 -E '^Sample .* stage-1 probe: ' "$log" | sed -e 's/^Sample //' -e 's/ ([^()]*_stage1_probe\*\.tsv)$//'
    done
fi
# --test_sparse: each level's wall and aligning time against the timed runs' medians, the reads given up, the records
# written against the timed runs', and the species against the last timed run's profile.
if [ $test_sparse = 1 ] && [ "$(wc -l < "$out/sparse.tsv")" -gt 1 ]; then
    echo
    awk -F'\t' "$median"'
        function med(x, k,   a, n, r) { split("", a); n = 0
            for (r = 1; r <= c[k]; r++) if (x[k, r] != "NA") a[++n] = x[k, r] + 0
            return n ? median(a, n) : 0 }
        function change(v, m) { return (m > 0 && v != "NA") ? sprintf("%+.1f%%", 100 * (v / m - 1)) : "NA" }
        function field(name) { return (name in col) ? $(col[name]) : "NA" }
        FNR == 1 { t = FILENAME ~ /sparse\.tsv$/ ? "sparse" : "runs"; if (t == "runs") for (i = 1; i <= NF; i++) col[$i] = i; next }
        t == "runs" { k = $1; r = ++c[k]; wall[k, r] = $5; aligning[k, r] = $11
                      cpu[k, r] = field("aligning_cpu_s"); reader[k, r] = field("reader_s_per_thread"); next }
        { k = $1; mw = med(wall, k); ma = med(aligning, k); mc = med(cpu, k); mr = med(reader, k)
          # (the comparisons before printf: a ">" among its arguments would redirect its output)
          dw = change($5, mw); da = change($9, ma); dc = change($29, mc); share = $15 == "NA" ? "NA" : sprintf("%.1f%%", 100 * $15)
          records = ($18 != "NA" && $19 != "NA" && $19 > 0) ? sprintf("%s of %s (%+.2f%%)", $18, $19, 100 * ($18 / $19 - 1)) : $18 " of " $19
          timed_reader = mr > 0 ? sprintf("%.2f", mr) : "NA"
          peak = $8 == "NA" ? "" : " " $8 "Gb"
          printf "%s --sparse_first %s: wall %s s%s (%s), aligning %s s (%s), given up %s, records %s, species %s (lost %s, gained %s), Bray-Curtis %s (species alone %s, unknown %s -> %s); aligning CPU %s s (%s) on %s cores, reader %s s a thread (timed runs: %s), gate %s s a thread\n",
              k, $4, $5, peak, dw, $9, da, share, records, $21, $23, $24, $25, $26, $27, $28, $29, dc, $30, $31, timed_reader, $32 }' \
        "$out/runs.tsv" "$out/sparse.tsv"
fi
[ -n "$grid_dropped" ] && { echo; echo "Note: $grid_dropped"; }
echo; echo "Results in $out (runs.tsv, cohort.tsv, stages.tsv,$([ -n "$sylph" ] && echo " sylph.tsv, sylph/,")$([ $stage1_probe = 1 ] && echo " probe.tsv, stage1_probe/,")$([ $test_sparse = 1 ] && echo " sparse.tsv, sparse_stages.tsv, sparse/,") environment.txt, logs/)."
