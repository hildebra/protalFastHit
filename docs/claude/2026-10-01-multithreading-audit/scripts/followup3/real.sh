#!/usr/bin/env bash
# Records straight into the .sam.zst (samroom: 88f20d4 + samroom.patch) against 88f20d4 (the samstep build, which
# only adds SAMSTEP lines on stderr), whole runs, niced:
#  1. -t 1 (records in a fixed order): paired-end 500k and 1k pairs, Nanopore 3 Mb, aligned and profiled by
#     each build: the SAMs' text (zstdcat) and every output but the SAM and the timings must be the same; zstd -t;
#     the SAMs' sizes (ls) and disk use (du).
#  2. -t 6, the 5M-pair sample, --no_profile, alternated, twice each: the SAM step's time and the run's.
set -uo pipefail
W=$HOME/mt-work/samreal; rm -rf $W; mkdir -p $W
declare -A BIN=( [base]=$HOME/mt-work/samstep/src/build/protal [room]=$HOME/mt-work/samroom/src/build/protal )
DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
declare -A ARGS=(
  [pe500k]="-1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz -2 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz"
  [pe1k]="-1 $P/rl150_p1000/sim/reads/rl150_p1000_s_1_R1.fq.gz -2 $P/rl150_p1000/sim/reads/rl150_p1000_s_1_R2.fq.gz"
  [ont3M]="-1 $P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz --read_type ont" )
for s in pe500k pe1k ont3M; do
  for b in base room; do
    nice -n 5 ${BIN[$b]} --db $DB ${ARGS[$s]} --prefix $s -o $W/$s.$b -t 1 --no_qcmsa > $W/$s.$b.log 2>&1 || echo "FAIL $s $b exit $?"
  done
  a=$W/$s.base/$s.sam.zst; c=$W/$s.room/$s.sam.zst
  cmp -s <(zstdcat $a) <(zstdcat $c) && echo "SAME SAM text $s" || echo "DIFFER SAM text $s"
  diff -r -q -x '*.sam.zst' -x '*_runtime.tsv' $W/$s.base $W/$s.room > $W/$s.diff && echo "SAME outputs $s ($(find $W/$s.room -type f | wc -l) files)" || { echo "DIFFER outputs $s"; head -5 $W/$s.diff; }
  zstd -tq $c && echo "zstd -t OK $s" || echo "FAIL zstd -t $s"
  echo "$s sizes: base $(stat -c %s $a) bytes, room $(stat -c %s $c) bytes, disk $(du -k $a | cut -f1) / $(du -k $c | cut -f1) KB; step $(grep -h 'Writing the SAM header' $W/$s.base.log | sed 's/.*took //') / $(grep -h 'Writing the SAM header' $W/$s.room.log | sed 's/.*took //')"
done
R=$HOME/bench071/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1
cat ${R}_R1.fq.gz ${R}_R2.fq.gz > /dev/null
for rep in 1 2; do for b in base room; do
  d=$W/deep; rm -rf $d
  nice -n 5 /usr/bin/time -f "%e" -o $W/time ${BIN[$b]} --db $DB -1 ${R}_R1.fq.gz -2 ${R}_R2.fq.gz --prefix s -o $d -t 6 --no_profile > $W/deep.log 2>&1
  echo "deep $b rep $rep: $(grep 'Aligning reads took' $W/deep.log | sed 's/.*took //'); SAM step $(grep 'Writing the SAM header' $W/deep.log | sed 's/.*took //'); wall $(cat $W/time) s; $(stat -c %s $d/s.sam.zst) bytes"
done; done
echo REAL DONE
