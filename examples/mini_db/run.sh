#!/usr/bin/env bash
# run.sh - reproducible end-to-end test of protal on a mini database.
#
#   bash examples/mini_db/run.sh [--rebuild] [WORKDIR]   (default WORKDIR: examples/mini_db/work)
#
# a) Build: simulate a sparse GTDB release (3 species x 3 genomes, seed 42), convert it
#    into a protal DB and run protal --build          -> WORKDIR/gtdb_r226, WORKDIR/protal_db
#    The DB is reused while the protal binary and the DB generator are unchanged;
#    --rebuild forces a new build.
# b) Test: simulate 30,000 read pairs from community.tsv (seed 7), profile them with
#    protal and compare with the truth (check_results.py)
#                                                     -> WORKDIR/reads, WORKDIR/run, WORKDIR/check.tsv
#
# All inputs are seeded: reruns regenerate identical files (WORKDIR/checksums.md5).
# WORKDIR/run_info.txt records the protal binary, git commit and settings.
#
# Environment: PROTAL  protal binary (default: build/protal, else protal on $PATH)
#              PYTHON  python 3 interpreter (default: python3; standard library only)
#              THREADS protal threads (default: 4)
# Needs ~3.5 GB of disk in WORKDIR: index.prx has a fixed-size key map.
# Exits 0 if all checks pass, 1 otherwise.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
tools="$root/scripts/mini_db"

# Options for simulate_gtdb_release.py (e.g. --lineages FILE); they are part of the DB stamp.
BUILD_ARGS=(--seed 42)
READS_SEED=7
READ_PAIRS=30000

rebuild=0
if [ "${1:-}" = "--rebuild" ]; then rebuild=1; shift; fi
work=${1:-$here/work}
python=${PYTHON:-python3}
threads=${THREADS:-4}

protal=${PROTAL:-}
if [ -z "$protal" ]; then
    if [ -x "$root/build/protal" ]; then protal="$root/build/protal"; else protal=$(command -v protal || true); fi
fi
if [ -z "$protal" ] || [ ! -x "$protal" ]; then
    echo "protal binary not found: build it (just baseline) or set PROTAL" >&2
    exit 1
fi
protal="$(cd "$(dirname "$protal")" && pwd)/$(basename "$protal")"

mkdir -p "$work"
work=$(cd "$work" && pwd)
db="$work/protal_db"
log() { printf '[mini_db] %s\n' "$*"; }
fail() { log "$1 failed, last lines of $2:"; tail -20 "$2" >&2; exit 1; }

# ---- a) build the database ------------------------------------------------------------------
stamp=$({ md5sum "$protal" "$tools"/{simulate_gtdb_release.py,gtdb_to_protal_db.py,markers_r226.tsv,build_mini_db.sh} \
          | cut -d' ' -f1; echo "build_args=${BUILD_ARGS[*]}"; } | md5sum | cut -d' ' -f1)
if [ "$rebuild" = 1 ] || [ ! -s "$db/index.prx" ] || [ "$(cat "$work/db.stamp" 2>/dev/null)" != "$stamp" ]; then
    log "a) building the mini DB in $work"
    rm -f "$work/db.stamp"
    PROTAL="$protal" PYTHON="$python" bash "$tools/build_mini_db.sh" "$work" "${BUILD_ARGS[@]}" \
        > "$work/build_mini_db.log" 2>&1 || fail "building the DB" "$work/build_mini_db.log"
    echo "$stamp" > "$work/db.stamp"
else
    log "a) reusing $db (same protal binary and generator; --rebuild to force)"
fi

# ---- b) simulate reads, profile them, compare with the truth ----------------------------------
reads="$work/reads"
run="$work/run"
rm -rf "$reads" "$run"
log "b) simulating $READ_PAIRS read pairs from $here/community.tsv"
"$python" "$tools/simulate_reads.py" --genomes "$work/gtdb_r226/simulation/genomes.tsv" \
    --community "$here/community.tsv" --out_prefix "$reads/mini" --pairs "$READ_PAIRS" --seed "$READS_SEED" \
    2> "$work/simulate_reads.log" || fail "simulating reads" "$work/simulate_reads.log"

log "   profiling with $protal"
"$protal" --db "$db" -1 "$reads/mini_R1.fq" -2 "$reads/mini_R2.fq" --prefix mini -o "$run" \
    -t "$threads" --no_qcmsa > "$work/protal.log" 2>&1 || fail "protal" "$work/protal.log"

{
    echo "date	$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "git_commit	$(git -C "$root" describe --always --dirty 2>/dev/null || echo unknown)"
    echo "protal	$protal"
    echo "protal_version	$("$protal" --version 2>/dev/null | tail -1)"
    echo "protal_md5	$(md5sum "$protal" | cut -d' ' -f1)"
    echo "python	$("$python" --version 2>&1)"
    echo "build_args	${BUILD_ARGS[*]}"
    echo "reads_seed	$READS_SEED"
    echo "read_pairs	$READ_PAIRS"
    echo "threads	$threads"
} > "$work/run_info.txt"
(cd "$work" && md5sum protal_db/{reference.fna,reference.map,internal_taxonomy.dmp,full_reference.fna,unique_kmers.tsv} \
    reads/{mini_R1.fq,mini_R2.fq,mini.truth.tsv}) > "$work/checksums.md5"

log "   checking the results"
"$python" "$here/check_results.py" --truth "$reads/mini.truth.tsv" --profile "$run/mini.profile" \
    --sam "$run/mini.sam" --taxonomy "$db/internal_taxonomy.dmp" --out "$work/check.tsv"
