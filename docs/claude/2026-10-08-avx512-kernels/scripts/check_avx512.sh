#!/usr/bin/env bash
# check_avx512.sh - the AVX-512 kernels (the flex scan of seed lookups, the syncmer scan) checked and timed on a node
# that has them (AMD Zen 4 or Intel Ice Lake and later), from a protal checkout built in build/ as Release with the
# unit tests (cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON; cmake --build build -j).
#
#   bash docs/claude/2026-10-08-avx512-kernels/scripts/check_avx512.sh OUT_DIR DB TYPE:R1[:R2] [TYPE:R1[:R2] ...]
#
# 1. The CPU's AVX-512 flags and the level protal's kernels use here (from the unit test that prints it).
# 2. The kernels' unit tests, which compare every level the CPU has with the scalar code: FlexScan.*, Syncmers.*,
#    PackedIndex.LookupsGiveTheSameSeedsAtEveryVectorLevel. A test that says "no AVX-512 here" did not test it.
# 3. Each sample once with PROTAL_SIMD=avx2 and once with the default (AVX-512 here): the SAM records (sorted) and
#    every other output file must be the same.
# 4. scripts/measure_performance.sh for both levels, alternated: OUT_DIR/perf_avx2_<i>, OUT_DIR/perf_avx512_<i>
#    (ROUNDS rounds, default 2, each REPEATS runs, default 2; COHORT=0). Compare the seeding and "taking the
#    k-mers" stages in their stages.tsv; the counts in runs.tsv must be equal.
#
# Environment: THREADS (default SLURM_CPUS_PER_TASK, else nproc), ROUNDS, REPEATS, PROTAL (default build/protal),
#              TESTS (default build/tests/protal_tests).
# On SLURM, a whole node, e.g.
#   sbatch --exclusive -c 32 --mem 100G --wrap "bash docs/claude/2026-10-08-avx512-kernels/scripts/check_avx512.sh \
#       avx512_r226 r226_db pe:R1.fq.gz:R2.fq.gz pb:hifi.fq.gz"
set -uo pipefail

[ $# -ge 3 ] || { sed -n '2,24p' "$0" >&2; exit 1; }
out=$1; db=$2; shift 2
threads=${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}
rounds=${ROUNDS:-2}
protal=${PROTAL:-build/protal}
tests=${TESTS:-build/tests/protal_tests}
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../.." && pwd)
mkdir -p "$out"
fail=0

echo "== 1. CPU ($(hostname), $(date -Is))"
grep -m1 'model name' /proc/cpuinfo
echo "AVX-512 flags: $(grep -m1 '^flags' /proc/cpuinfo | tr ' ' '\n' | grep -E '^avx512' | sort | tr '\n' ' ')"
for f in avx512f avx512bw avx512vl avx512dq avx512vbmi avx512_vbmi2 avx512_vpopcntdq; do
    grep -m1 '^flags' /proc/cpuinfo | tr ' ' '\n' | grep -qx "$f" || echo "missing: $f (protal's kernels use AVX2 here)"
done

echo "== 2. unit tests of the kernels"
"$tests" --gtest_filter='FlexScan.*:Syncmers.*:PackedIndex.LookupsGiveTheSameSeedsAtEveryVectorLevel' > "$out/unit.log" 2>&1
code=$?
grep -E "this CPU:|no AVX|^\[  (PASSED|FAILED)|FAILED  \]|blocks, " "$out/unit.log"
[ $code = 0 ] || { echo "unit tests FAILED (exit $code): $out/unit.log"; fail=1; }
grep -q "this CPU: AVX-512" "$out/unit.log" || echo "NOTE: this CPU does not have protal's AVX-512 level; steps 3 and 4 compare AVX2 with AVX2"

echo "== 3. the same outputs at both levels"
for spec in "$@"; do
    IFS=: read -r type r1 r2 <<< "$spec"
    args=(-1 "$r1"); [ -n "$r2" ] && args+=(-2 "$r2")
    for level in avx2 auto; do
        o=$out/cmp/$type.$level; rm -rf "$o"; mkdir -p "$out/cmp"
        PROTAL_SIMD=$level "$protal" --db "$db" "${args[@]}" --read_type "$type" --prefix s -o "$o" -t "$threads" --no_qcmsa \
            > "$o.log" 2>&1 || { echo "$type $level: protal failed, $o.log"; fail=1; }
    done
    a=$out/cmp/$type.avx2; b=$out/cmp/$type.auto
    sam() { if [ -f "$1/s.sam.zst" ]; then zstd -dc "$1/s.sam.zst"; else zcat "$1/s.sam.gz"; fi | grep -v '^@PG\|^@CO' | LC_ALL=C sort | md5sum | cut -c1-16; }
    sa=$(sam "$a"); sb=$(sam "$b")
    diffs=$(cd "$a" && find . -type f ! -name '*.sam.zst' ! -name '*.sam.gz' ! -name '*.err' ! -name '*_runtime.tsv' |
            while read -r f; do cmp -s "$f" "$b/$f" || echo "$f"; done | tr '\n' ' ')
    if [ "$sa" = "$sb" ] && [ -z "$diffs" ]; then echo "$type: identical (SAM records and every other file)"
    else echo "$type: DIFFER: SAM $sa vs $sb; other files: ${diffs:-none}"; fail=1; fi
done

echo "== 4. timing, $rounds rounds alternated"
for i in $(seq 1 "$rounds"); do
    for level in avx2 avx512; do
        v=$level; [ $level = avx512 ] && v=auto
        PROTAL_SIMD=$v PROTAL=$protal THREADS=$threads REPEATS=${REPEATS:-2} COHORT=0 \
            bash "$repo/scripts/measure_performance.sh" "$out/perf_${level}_$i" "$db" "$@" > "$out/perf_${level}_$i.log" 2>&1
        echo "round $i $level: $(tail -3 "$out/perf_${level}_$i.log" | tr '\n' ' ')"
    done
done
echo "== done: $([ $fail = 0 ] && echo "all checks passed" || echo "SOME CHECKS FAILED"); outputs in $out"
exit $fail
