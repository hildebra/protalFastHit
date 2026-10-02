#!/usr/bin/env bash
# WFA2's AVX2 extend kernels chosen at run time (wfa: f07acb3 + wfa.patch) against f07acb3's code with the portable
# kernels (samroomhead: 5c7d64c + 8798733's change, the same src), whole runs, niced:
#  1. -t 1 (records in a fixed order): paired-end 500k and 1k pairs, Nanopore 3 Mb: the SAM text and every output
#     but the SAM and the timings must be the same.
#  2. -t 6, alternated, three times: the 5M-pair sample and Nanopore 90 Mb, aligned and profiled by the old build,
#     the new baseline build and the new protal_avx2; the stage times, and the records of each SAM sorted (their
#     order changes from run to run at 6 threads) must be the same in every run.
# Rows into ~/mt-work/wfareal/runs.tsv.
set -uo pipefail
W=$HOME/mt-work/wfareal; rm -rf $W; mkdir -p $W
declare -A BIN=( [old]=$HOME/mt-work/samroomhead/src/build/protal [new]=$HOME/mt-work/wfa/src/build/protal
                 [new_avx2]=$HOME/mt-work/wfa/src/build/protal_avx2 )
DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points; D=$HOME/bench071/samples_deep/points
declare -A ARGS=(
  [pe500k]="-1 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R1.fq.gz -2 $P/rl150_p500000/sim/reads/rl150_p500000_s_1_R2.fq.gz"
  [pe1k]="-1 $P/rl150_p1000/sim/reads/rl150_p1000_s_1_R1.fq.gz -2 $P/rl150_p1000/sim/reads/rl150_p1000_s_1_R2.fq.gz"
  [ont3M]="-1 $P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz --read_type ont"
  [pe5M]="-1 $D/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz -2 $D/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R2.fq.gz"
  [ont90M]="-1 $P/ont_b90000000/sim/reads/ont_b90000000_s_1.fq.gz --read_type ont" )
for s in pe500k pe1k ont3M; do
  for b in old new; do
    nice -n 5 ${BIN[$b]} --db $DB ${ARGS[$s]} --prefix $s -o $W/$s.$b -t 1 --no_qcmsa > $W/$s.$b.log 2>&1 || echo "FAIL $s $b exit $?"
  done
  cmp -s <(zstdcat $W/$s.old/$s.sam.zst) <(zstdcat $W/$s.new/$s.sam.zst) && echo "SAME SAM text $s" || echo "DIFFER SAM text $s"
  diff -r -q -x '*.sam.zst' -x '*_runtime.tsv' $W/$s.old $W/$s.new > $W/$s.diff && echo "SAME outputs $s ($(find $W/$s.new -type f | wc -l) files)" || { echo "DIFFER outputs $s"; head -5 $W/$s.diff; }
  rm -rf $W/$s.old $W/$s.new
done
seconds() { awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'; }
printf "sample\tbinary\trep\taligning_s\tprofiling_s\twall_s\tuser_s\tload\trecords_md5\n" > $W/runs.tsv
for rep in 1 2 3; do
  for s in pe5M ont90M; do
    for b in old new new_avx2; do
      d=$W/run; rm -rf $d
      load=$(cut -d' ' -f1 /proc/loadavg)
      nice -n 5 /usr/bin/time -f "%e %U" -o $W/time ${BIN[$b]} --db $DB ${ARGS[$s]} --prefix $s -o $d -t 6 --no_qcmsa > $W/log 2>&1 || echo "FAIL $s $b exit $?"
      md5=$(zstdcat $d/$s.sam.zst | grep -v '^@' | LC_ALL=C sort | md5sum | cut -c1-12)
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $s $b $rep $(grep 'Aligning reads took' $W/log | sed 's/.*took //' | seconds) \
          $(grep 'Profiling took' $W/log | sed 's/.*took //' | seconds) "$(tail -1 $W/time | tr ' ' '\t')" $load $md5 | tee -a $W/runs.tsv
    done
  done
done
echo WFAREAL DONE
