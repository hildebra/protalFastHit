#!/usr/bin/env bash
# Compares two output folders of runs.sh: the SAMs as sorted sets of records, the other outputs with diff -r (left out:
# the timings *_runtime.tsv and the per-taxon *.statistics.tsv). Prints "SAM same|DIFF (n records)  outputs same|DIFF".
#   compare.sh <dir a> <dir b> [label]
set -uo pipefail
A=$1; B=$2; L=${3:-$A vs $B}
sam() { local f; f=$(ls $1/*.sam.zst $1/*.sam 2>/dev/null | head -1); case $f in *.zst) zstd -dcq $f;; *) cat $f;; esac | grep -v '^@' | LC_ALL=C sort; }
n=$(diff <(sam $A) <(sam $B) | grep -c '^[<>]')
s=$([ "$n" = 0 ] && echo "SAM same" || echo "SAM DIFF ($n records)")
d=$(diff -r -x '*_runtime.tsv' -x '*.sam.zst' -x '*.sam' -x '*.statistics.tsv' $A $B 2>&1)
o=$([ -z "$d" ] && echo "outputs same" || echo "outputs DIFF ($(echo "$d" | grep -c '^[<>]') lines in $(echo "$d" | grep -c '^diff') files)")
echo "$L: $s, $o"
