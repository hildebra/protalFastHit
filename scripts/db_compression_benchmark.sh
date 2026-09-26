#!/usr/bin/env bash
# db_compression_benchmark.sh - how well does zstd compress a protal database, and how fast does
# protal load it, raw and compressed, from where the database lives (e.g. network storage)?
#
#   bash scripts/db_compression_benchmark.sh compress DB_DIR WORK_DIR
#   bash scripts/db_compression_benchmark.sh load DB_DIR [DB_DIR ...]
#
# compress  For index.prx and reference.fna of DB_DIR (raw or .zst), and every LEVELS x WINDOW_LOGS
#           setting: compressed size, ratio, compression and decompression time. Also the time to
#           read the raw file. WORK_DIR (ideally on the same storage as DB_DIR) receives the
#           compressed files and WORK_DIR/db_zst: a copy of the database with the best-ratio
#           setting (other files symlinked), ready for the load mode.
# load      Runs protal on a few reads against each DB_DIR and reports "Load Index" and
#           "Preload genomes" times. The page cache serves files read recently; for cold-read
#           numbers run each database once on a node that has not read it (e.g. a fresh job).
#
# Environment: LEVELS (default "12 19"), WINDOW_LOGS (default "27"; 0 = no long-distance matching),
#              THREADS (default: nproc), PROTAL (default: build/protal, load mode).
# Needs the zstd CLI. Compressing a ~50 GB database at level 19 takes roughly
# 50 GB / (3 MB/s x THREADS), e.g. ~20 min with 16 threads.
set -euo pipefail

mode=${1:-}
levels=${LEVELS:-12 19}
windows=${WINDOW_LOGS:-27}
threads=${THREADS:-$(nproc)}
protal=${PROTAL:-build/protal}

usage() { sed -n '2,9p' "$0" >&2; exit 1; }
command -v zstd > /dev/null || { echo "the zstd CLI is needed" >&2; exit 1; }
now() { date +%s.%N; }
elapsed() { awk -v a="$1" -v b="$2" 'BEGIN { printf "%.1f", b - a }'; }
mb() { awk -v b="$1" 'BEGIN { printf "%.1f", b / 1048576 }'; }

compress_mode() {
    local db=$1 work=$2
    mkdir -p "$work"
    printf '%-14s %5s %6s %12s %12s %7s %10s %10s %10s\n' \
        file level window raw_MB zst_MB ratio comp_s decomp_s read_raw_s
    local best_db="$work/db_zst"
    rm -rf "$best_db" && mkdir -p "$best_db"
    for f in "$db"/*; do
        case "$(basename "$f")" in index.prx*|reference.fna*) ;; *) ln -s "$(cd "$(dirname "$f")" && pwd)/$(basename "$f")" "$best_db/";; esac
    done
    for name in index.prx reference.fna; do
        local raw="$db/$name"
        if [ ! -e "$raw" ]; then
            [ -e "$raw.zst" ] || { echo "neither $raw nor $raw.zst exists" >&2; exit 1; }
            echo "decompressing $raw.zst to $work/$name for the benchmark" >&2
            zstd -q -d -f --long=31 -T"$threads" "$raw.zst" -o "$work/$name"
            raw="$work/$name"
        fi
        local size; size=$(stat -Lc%s "$raw")
        local t0; t0=$(now); cat "$raw" > /dev/null; local read_s; read_s=$(elapsed "$t0" "$(now)")
        local best="" best_size=0
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
                if [ -z "$best" ] || [ "$zsize" -lt "$best_size" ]; then best=$out; best_size=$zsize; fi
            done
        done
        ln -s "$(cd "$(dirname "$best")" && pwd)/$(basename "$best")" "$best_db/$name.zst"
    done
    echo
    echo "Best-ratio compressed copy: $best_db (compare with: $0 load $db $best_db)"
}

load_mode() {
    [ -x "$protal" ] || { echo "protal binary not found at $protal (set PROTAL)" >&2; exit 1; }
    local tmp; tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' RETURN
    # 20 read pairs from the first reference gene of at least 300 bp.
    local ref="$1/reference.fna"
    local gene
    if [ -e "$ref" ]; then gene=$(awk '!/^>/ && length($0) >= 300 { print; exit }' "$ref")
    else gene=$( (zstd -dc "$ref.zst" || true) | awk '!/^>/ && length($0) >= 300 { print; exit }'); fi
    [ -n "$gene" ] || { echo "no reference gene of >= 300 bp found in $ref" >&2; exit 1; }
    awk -v g="$gene" 'BEGIN {
        for (i = 0; i < 20; i++) {
            s = substr(g, 1 + (i * 7) % (length(g) - 299), 300); q = sprintf("%100s", ""); gsub(/ /, "I", q)
            r2 = ""; t = substr(s, 201, 100)
            for (j = 100; j >= 1; j--) { c = substr(t, j, 1); r2 = r2 (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : "A") }
            print "@r" i "\n" substr(s, 1, 100) "\n+\n" q > "'"$tmp"'/r_R1.fq"
            print "@r" i "\n" r2 "\n+\n" q > "'"$tmp"'/r_R2.fq"
        } }'
    printf '%-60s %16s %16s\n' database load_index preload_genomes
    local n=0
    for db in "$@"; do
        n=$((n + 1))
        local log="$tmp/run$n.log"
        "$protal" --db "$db" -1 "$tmp/r_R1.fq" -2 "$tmp/r_R2.fq" --prefix r -o "$tmp/out$n" \
            -t "$threads" --no_qcmsa --no_strains > "$log" 2>&1 || { echo "protal failed on $db:" >&2; tail -5 "$log" >&2; continue; }
        printf '%-60s %16s %16s\n' "$db" "$(grep -m1 'Load Index took' "$log" | sed 's/.*took //')" \
            "$(grep -m1 'Preload genomes took' "$log" | sed 's/.*took //')"
    done
}

case "$mode" in
    compress) [ $# -eq 3 ] || usage; compress_mode "$2" "$3" ;;
    load) [ $# -ge 2 ] || usage; shift; load_mode "$@" ;;
    *) usage ;;
esac
