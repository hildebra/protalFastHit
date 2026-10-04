#!/usr/bin/env bash
# head against the working tree (and the working tree without the k-mer screen): outputs identical? times, counts.
set -uo pipefail
W=$HOME/perf-gtdb; H=$W/head/build/protal; N=$W/work/build/protal; DB=$HOME/bench071/V073/protal_db
P=$HOME/bench071/samples/points
R1=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz; R2=$P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz
D1=$HOME/bench071/samples_deep/points/rl150_p5000000_s_1_R1.fq.gz; D2=$HOME/bench071/samples_deep/points/rl150_p5000000_s_1_R2.fq.gz
PB=$(ls $HOME/bench071/samples_lr073/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz 2>/dev/null || ls $HOME/bench071/samples/points/pb_b90000000/sim/reads/pb_b90000000_s_1.fq.gz)
echo "pb reads: $PB"
run() { # name binary threads outdir reads...
  local name=$1 bin=$2 t=$3 o=$4; shift 4
  rm -rf $o
  /usr/bin/time -f "%e s wall, %U s user, %M KB" $bin --db $DB "$@" --prefix s -o $o -t $t --no_qcmsa > $o.log 2> $o.time
  echo "$name: $(tail -1 $o.time); $(grep -E '^Aligning reads took' $o.log); $(grep -E '^Profiling took' $o.log)"
  grep -E '^Sample s: .*candidate|^Profiling sample' $o.log | sed 's/^/    /'
}
same() { # a b
  if cmp -s <(zstd -dc $1/s.sam.zst) <(zstd -dc $2/s.sam.zst); then echo "    SAM identical ($(zstd -dc $2/s.sam.zst | grep -vc '^@') records)"; else echo "    SAM DIFFERS"; diff <(zstd -dc $1/s.sam.zst) <(zstd -dc $2/s.sam.zst) | head -6; fi
  if diff -r -x '*_runtime.tsv' -x '*.sam.zst' $1 $2 > $W/diff.txt; then echo "    other outputs identical"; else echo "    OUTPUTS DIFFER"; head -8 $W/diff.txt; fi
}
echo "== 500k pairs, 1 thread"
run head $H 1 $W/s.head -1 $R1 -2 $R2 --read_type pe
run work $N 1 $W/s.work -1 $R1 -2 $R2 --read_type pe
run work-noscreen $N 1 $W/s.work0 -1 $R1 -2 $R2 --read_type pe --no_alignment_screen
same $W/s.head $W/s.work; same $W/s.head $W/s.work0
echo "== PacBio 90 Mb, 6 threads"
run head $H 6 $W/p.head -1 $PB --read_type pb
run work $N 6 $W/p.work -1 $PB --read_type pb
run work-noscreen $N 6 $W/p.work0 -1 $PB --read_type pb --no_alignment_screen
same $W/p.head $W/p.work; same $W/p.head $W/p.work0
echo "== 5M pairs, 6 threads (two rounds)"
for round in 1 2; do
  run "head $round" $H 6 $W/d.head.$round -1 $D1 -2 $D2 --read_type pe
  run "work $round" $N 6 $W/d.work.$round -1 $D1 -2 $D2 --read_type pe
done
same $W/d.head.1 $W/d.work.1
echo "== runtime table of the deep work run"; cat $W/d.work.1/misc/s_runtime.tsv
