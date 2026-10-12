#!/usr/bin/env bash
# Exact instruction counts of the paired-end loop with and without --sparse_first (docs/claude/2026-10-12-sparse-first):
# callgrind on one thread, collecting inside RunPairedEnd's parallel region only (its OpenMP body: --toggle-collect
# flips at every entry, so a pattern that also matched RunPairedEnd itself would switch it off again inside), aligning
# only; the loop's total and its top functions.
#   bash sf_callgrind.sh PROTAL DATABASE R1 R2 WORK ARM [protal options...]
set -uo pipefail
protal=$1 db=$2 r1=$3 r2=$4 work=$5 arm=$6; shift 6
mkdir -p "$work"
rm -rf "$work/run_$arm"
taskset -c 0-3 nice -n 5 valgrind --tool=callgrind --callgrind-out-file="$work/cg_$arm.out" \
    --toggle-collect='*RunPairedEnd*_omp_fn*' "$protal" --db "$db" -1 "$r1" -2 "$r2" --prefix s -o "$work/run_$arm" -t 1 \
    --no_strains --no_profile "$@" > "$work/cg_$arm.log" 2>&1
echo "$arm exit $?"
callgrind_annotate --inclusive=yes "$work/cg_$arm.out" 2> /dev/null | head -70 > "$work/cg_$arm.txt"
callgrind_annotate "$work/cg_$arm.out" 2> /dev/null | head -45 > "$work/cg_${arm}_self.txt"
grep -m1 "PROGRAM TOTALS" "$work/cg_$arm.txt"
sed -n '/file:function/,$p' "$work/cg_${arm}_self.txt" | head -30 | cut -c1-200
