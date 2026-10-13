#!/usr/bin/env bash
# fi_check_avx512.sh - --flex_indels' AVX-512 kernel (flex_scan::ScanAvx512<true>) checked on a node that has protal's
# AVX-512 level (AMD Zen 4, Intel Ice Lake and later; the laptop it was written on has none), from a protal checkout
# built in build/ as Release with the unit tests (cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON).
#
#   bash docs/claude/2026-10-13-flex-indels/scripts/fi_check_avx512.sh OUT_DIR DB TYPE:R1[:R2] [TYPE:R1[:R2] ...]
#
# 1. The unit tests of the scan and the lookups, which compare every level the CPU has with the scalar code, with and
#    without --flex_indels' scores (FlexScan.*, PackedIndex.*). "this CPU: AVX-512" must be printed, else nothing was
#    tested.
# 2. Each sample with --flex_indels all, once with PROTAL_SIMD=avx2 and once with the default (AVX-512 here): the SAM
#    records (sorted) and the profile must be the same. Then the default level without the option, for the seeding
#    stage's time with and without it (misc/s_runtime.tsv's "Seeding" row).
# Environment: THREADS (default SLURM_CPUS_PER_TASK, else nproc), PROTAL (default build/protal), TESTS (default
# build/tests/protal_tests). On SLURM a whole node, e.g.
#   sbatch --exclusive -c 32 --mem 100G --wrap "bash docs/claude/2026-10-13-flex-indels/scripts/fi_check_avx512.sh \
#       flex_r226 r226_db ont:soil_ont.fq.gz pe:R1.fq.gz:R2.fq.gz"
set -uo pipefail
[ $# -ge 3 ] || { sed -n '2,19p' "$0" >&2; exit 1; }
out=$1; db=$2; shift 2
threads=${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}
protal=${PROTAL:-build/protal}
tests=${TESTS:-build/tests/protal_tests}
mkdir -p "$out"
fail=0

echo "== 1. unit tests ($(hostname), $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2))"
"$tests" --gtest_filter='FlexScan.*:PackedIndex.*' > "$out/unit.log" 2>&1
code=$?
grep -E "this CPU:|no AVX|^\[  (PASSED|FAILED)|FAILED  \]|blocks, " "$out/unit.log"
[ $code = 0 ] || { echo "unit tests FAILED (exit $code): $out/unit.log"; fail=1; }
grep -q "this CPU: AVX-512" "$out/unit.log" || echo "NOTE: no AVX-512 level here: the AVX-512 kernel was not tested"

echo "== 2. --flex_indels all at both levels, and the seeding's time without and with it"
sam() { zstd -dc "$1/s.sam.zst" | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-16; }
seeding() { awk -F'\t' '$1 == "Seeding" { printf "%.1f s (%.2f a thread)", $2, $4 }' "$1/misc/s_runtime.tsv"; }
for spec in "$@"; do
    IFS=: read -r type r1 r2 <<< "$spec"
    args=(-1 "$r1"); [ -n "$r2" ] && args+=(-2 "$r2")
    for arm in avx2 auto plain; do
        o=$out/runs/$type.$arm; rm -rf "$o"; mkdir -p "$out/runs"
        level=$arm extra=(--flex_indels all); [ $arm = plain ] && { level=auto; extra=(); }
        PROTAL_SIMD=$level "$protal" --db "$db" "${args[@]}" --read_type "$type" --prefix s -o "$o" -t "$threads" --no_qcmsa \
            --no_strains "${extra[@]}" > "$o.log" 2>&1 || { echo "$type $arm: protal failed, $o.log"; fail=1; }
    done
    a=$out/runs/$type.avx2 b=$out/runs/$type.auto
    if [ "$(sam "$a")" = "$(sam "$b")" ] && cmp -s "$a/s.profile" "$b/s.profile"; then
        echo "$type: AVX2 and AVX-512 the same with --flex_indels (records $(sam "$a"))"
    else
        echo "$type: AVX2 and AVX-512 DIFFER with --flex_indels"; fail=1
    fi
    echo "$type seeding: without $(seeding "$out/runs/$type.plain"), with --flex_indels $(seeding "$b")"
done
[ $fail = 0 ] && echo "== all checks passed" || echo "== a check FAILED"
exit $fail
