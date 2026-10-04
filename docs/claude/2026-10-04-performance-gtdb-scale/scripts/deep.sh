#!/usr/bin/env bash
set -uo pipefail
W=$HOME/perf-gtdb; H=$W/head/build/protal; N=$W/work/build/protal; DB=$HOME/bench071/V073/protal_db
D1=$(find $HOME/bench071 -name 'rl150_p5000000_s_1_R1.fq.gz' | head -1); D2=${D1/_R1/_R2}
echo "deep reads: $D1"
run() { local name=$1 bin=$2 t=$3 o=$4; shift 4; rm -rf $o
  /usr/bin/time -f "%e s wall, %U s user, %M KB" $bin --db $DB "$@" --prefix s -o $o -t $t --no_qcmsa > $o.log 2> $o.time
  echo "$name: $(tail -1 $o.time); $(grep -E '^Aligning reads took' $o.log); $(grep -E '^Profiling took' $o.log)"
  grep -E '^Sample s: .*candidate|^Profiling sample' $o.log | sed 's/^/    /'; }
for round in 1 2; do
  run "head $round" $H 6 $W/d.head.$round -1 $D1 -2 $D2 --read_type pe
  run "work $round" $N 6 $W/d.work.$round -1 $D1 -2 $D2 --read_type pe
done
run "work-noscreen" $N 6 $W/d.work0 -1 $D1 -2 $D2 --read_type pe --no_alignment_screen
a=$W/d.head.1; b=$W/d.work.1
if cmp -s <(zstd -dc $a/s.sam.zst | sort) <(zstd -dc $b/s.sam.zst | sort); then echo "SAM identical as sets ($(zstd -dc $b/s.sam.zst | grep -vc '^@') records)"; else echo "SAM DIFFERS"; fi
if diff -r -x '*_runtime.tsv' -x '*.sam.zst' $a $b > $W/ddiff.txt; then echo "other outputs identical"; else echo "OUTPUTS DIFFER"; head -5 $W/ddiff.txt | cut -c1-200; fi
echo "== runtime table of the deep work run"; cat $b/misc/s_runtime.tsv
