#!/usr/bin/env bash
# db_compression_benchmark.sh - how well does zstd compress a protal database, and how fast does
# protal load it, raw and compressed, from where the database lives (e.g. network storage)?
#
#   bash scripts/db_compression_benchmark.sh compress DB_DIR WORK_DIR
#   bash scripts/db_compression_benchmark.sh load DB_DIR [DB_DIR ...]
#
# compress  For index.prx and reference.fna of DB_DIR, and every LEVELS x WINDOW_LOGS setting:
#           compressed size, ratio, compression and decompression time (zstd CLI, one frame, i.e.
#           the best ratio). Also the time to read the raw file. A compressed database (.zst files
#           or database.protal) is first written raw into WORK_DIR/db_raw (protal --decompress_db).
#           Then WORK_DIR/db_single: the database as protal writes it (database.protal, seekable
#           64 MB frames, loaded with -t threads; `protal --compress_db`, level SEEKABLE_LEVEL) with
#           its size and time. WORK_DIR should be on the same storage as DB_DIR.
# load      Runs protal on a few reads against each DB_DIR (a folder or a database.protal), with
#           each thread count of THREAD_LIST, and reports "Load Index" and "Preload genomes" times.
#           The page cache serves files read recently; for cold-read numbers run each database
#           once on a node that has not read it (e.g. a fresh job).
#
# Environment: LEVELS (default "12 19"), WINDOW_LOGS (default "27"; 0 = no long-distance matching),
#              THREADS (default: nproc), THREAD_LIST (load; default "1 THREADS"),
#              SEEKABLE_LEVEL (default 19), FRAME_MB (default 64), PROTAL (default: build/protal).
# Needs the zstd CLI. Compressing a ~50 GB database at level 19 takes roughly
# 50 GB / (3 MB/s x THREADS), e.g. ~20 min with 16 threads (per setting).
set -euo pipefail

mode=${1:-}
levels=${LEVELS:-12 19}
windows=${WINDOW_LOGS:-27}
threads=${THREADS:-$(nproc)}
thread_list=${THREAD_LIST:-1 $threads}
protal=${PROTAL:-build/protal}

usage() { sed -n '2,20p' "$0" >&2; exit 1; }
command -v zstd > /dev/null || { echo "the zstd CLI is needed" >&2; exit 1; }
now() { date +%s.%N; }
elapsed() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.1f", b - a }'; }
mb() { awk -v b="$1" 'BEGIN { printf "%.1f", b / 1048576 }'; }
abs() { echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"; }

compress_mode() {
    local db=$1 work=$2
    mkdir -p "$work"
    if [ ! -e "$db/index.prx" ] || [ ! -e "$db/reference.fna" ]; then
        # The column-format index.prx.zst and database.protal do not zstd -d to raw files.
        [ -x "$protal" ] || { echo "PROTAL is needed to write $db raw" >&2; exit 1; }
        local raw_db="$work/db_raw"
        rm -rf "$raw_db" && mkdir -p "$raw_db"
        for f in "$db"/*; do ln -s "$(abs "$f")" "$raw_db/"; done
        echo "writing $db raw into $raw_db for the benchmark (protal --decompress_db)" >&2
        "$protal" --decompress_db --db "$raw_db" -t "$threads" > "$work/decompress_db.log" 2>&1 \
            || { echo "protal --decompress_db failed:" >&2; tail -5 "$work/decompress_db.log" >&2; exit 1; }
        db=$raw_db
    fi
    printf '%-14s %5s %6s %12s %12s %7s %10s %10s %10s\n' \
        file level window raw_MB zst_MB ratio comp_s decomp_s read_raw_s
    for name in index.prx reference.fna; do
        local raw="$db/$name"
        local size; size=$(stat -Lc%s "$raw")
        local t0; t0=$(now); cat "$raw" > /dev/null; local read_s; read_s=$(elapsed "$t0" "$(now)")
        for level in $levels; do
            for w in $windows; do
                local out="$work/$name.L$level.W$w.zst" flags=(-q -f -T"$threads" "-$level")
                [ "$level" -gt 19 ] && flags+=(--ultra)
                [ "$w" -gt 0 ] && flags+=("--long=$w")
                t0=$(now); zstd "${flags[@]}" "$raw" -o "$out"; local comp_s; comp_s=$(elapsed "$t0" "$(now)")
                t0=$(now); zstd -q -dc --long=31 "$out" > /dev/null; local decomp_s; decomp_s=$(elapsed "$t0" "$(now)")
                local zsize; zsize=$(stat -c%s "$out")
                printf '%-14s %5s %6s %12s %12s %7s %10s %10s %10s\n' "$name" "$level" "$w" "$(mb "$size")" \
                    "$(mb "$zsize")" "$(awk -v a="$size" -v b="$zsize" 'BEGIN { printf "%.2f", a / b }')" \
                    "$comp_s" "$decomp_s" "$read_s"
            done
        done
    done

    [ -x "$protal" ] || { echo; echo "PROTAL not found ($protal): no single-file copy"; return; }
    local single="$work/db_single"
    rm -rf "$single" && mkdir -p "$single"
    for f in "$db"/*; do ln -s "$(abs "$f")" "$single/"; done
    echo
    local t0; t0=$(now)
    "$protal" --compress_db --db "$single" -t "$threads" --compress_level "${SEEKABLE_LEVEL:-19}" \
        --compress_frame_mb "${FRAME_MB:-64}" > "$work/compress_db.log" 2>&1 \
        || { echo "protal --compress_db failed:"; tail -5 "$work/compress_db.log"; return; }
    echo "protal --compress_db (level ${SEEKABLE_LEVEL:-19}, ${FRAME_MB:-64} MB frames, $threads threads): $(elapsed "$t0" "$(now)") s"
    grep -E "(\.zst|database\.protal): " "$work/compress_db.log" || true
    echo "Compare load times with: $0 load $db $single"
}

# A reference gene of >= 300 bp from the first database given that has one to read: a separate
# reference.fna(.zst), else a database.protal (its parts decompressed in turn until a gene shows up).
reference_gene() {
    local db ref
    for db in "$@"; do
        ref="$db/reference.fna"
        if [ -e "$ref" ]; then awk '!/^>/ && length($0) >= 300 { print; exit }' "$ref"; return; fi
        if [ -e "$ref.zst" ]; then (zstd -dc "$ref.zst" || true) | awk '!/^>/ && length($0) >= 300 { print; exit }'; return; fi
    done
    for db in "$@"; do
        [ -f "$db" ] && ref=$db || ref="$db/database.protal"
        if [ -e "$ref" ]; then
            (zstd -dc "$ref" 2> /dev/null || true) | LC_ALL=C awk 'length($0) >= 300 && /^[ACGTN]+$/ { print; exit }'
            return
        fi
    done
}

load_mode() {
    [ -x "$protal" ] || { echo "protal binary not found at $protal (set PROTAL)" >&2; exit 1; }
    local tmp; tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' RETURN
    # 20 read pairs from a reference gene of at least 300 bp.
    local gene
    gene=$(reference_gene "$@")
    [ -n "$gene" ] || { echo "no reference gene of >= 300 bp found in $*" >&2; exit 1; }
    awk -v g="$gene" 'BEGIN {
        for (i = 0; i < 20; i++) {
            s = substr(g, 1 + (i * 7) % (length(g) - 299), 300); q = sprintf("%100s", ""); gsub(/ /, "I", q)
            r2 = ""; t = substr(s, 201, 100)
            for (j = 100; j >= 1; j--) { c = substr(t, j, 1); r2 = r2 (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : "A") }
            print "@r" i "\n" substr(s, 1, 100) "\n+\n" q > "'"$tmp"'/r_R1.fq"
            print "@r" i "\n" r2 "\n+\n" q > "'"$tmp"'/r_R2.fq"
        } }'
    printf '%-60s %8s %16s %16s\n' database threads load_index preload_genomes
    local n=0
    for db in "$@"; do
        for t in $thread_list; do
            n=$((n + 1))
            local log="$tmp/run$n.log"
            "$protal" --db "$db" -1 "$tmp/r_R1.fq" -2 "$tmp/r_R2.fq" --prefix r -o "$tmp/out$n" \
                -t "$t" --no_qcmsa --no_strains > "$log" 2>&1 || { echo "protal failed on $db:" >&2; tail -5 "$log" >&2; continue; }
            printf '%-60s %8s %16s %16s\n' "$db" "$t" "$(grep -m1 'Load Index took' "$log" | sed 's/.*took //')" \
                "$(grep -m1 'Preload genomes took' "$log" | sed 's/.*took //')"
        done
    done
}

case "$mode" in
    compress) [ $# -eq 3 ] || usage; compress_mode "$2" "$3" ;;
    load) [ $# -ge 2 ] || usage; shift; load_mode "$@" ;;
    *) usage ;;
esac
