#!/usr/bin/env bash
# callgrind (instructions, one thread) of four builds of d11381f's code on isa_cg.sh's three workloads:
# base = plain x86-64 (oneb-ref/protal), v3 = protal_avx2, cl = the target_clones build (oneb/protal),
# v2 = x86-64-v2 baseline (oneb-v2/protal). Four at a time, niced. Per-function self instructions in
# ~/mt-work/onebcg/<workload>.<build>.txt, then whether cl, v2 and v3 give base's outputs.
set -uo pipefail
W=$HOME/mt-work/onebcg; mkdir -p $W
I=$HOME/mt-work/isacg
DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
declare -A BIN=( [base]=$HOME/mt-work/oneb-ref/src/build/protal [v3]=$HOME/mt-work/oneb-ref/src/build/protal_avx2
                 [cl]=$HOME/mt-work/oneb/src/build/protal [v2]=$HOME/mt-work/oneb-v2/src/build/protal )
declare -A ARGS=(
  [pe100k]="-1 $I/r1.fq -2 $I/r2.fq --no_profile"
  [ont3M]="-1 $P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz --read_type ont --no_profile"
  [prof500k]="--profile_only $HOME/bench071/runs/v071.full.pe.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst" )
one() {  # one WORKLOAD BUILD
  local w=$1 b=$2
  rm -rf $W/out.$w.$b
  nice valgrind --tool=callgrind --callgrind-out-file=$W/cg.$w.$b ${BIN[$b]} --db $DB ${ARGS[$w]} --prefix s -o $W/out.$w.$b -t 1 --no_qcmsa > $W/log.$w.$b 2>&1
  callgrind_annotate --threshold=100 $W/cg.$w.$b 2>/dev/null | c++filt > $W/$w.$b.txt
  echo "done $w $b: $(grep 'PROGRAM TOTALS' $W/$w.$b.txt)"
}
for w in pe100k ont3M prof500k; do
  for b in base v3 cl v2; do one $w $b & done
  wait
done
for w in pe100k ont3M prof500k; do
  for b in cl v2 v3; do
    r=$(diff -r -q -x '*_runtime.tsv' -x '*.sam.zst' $W/out.$w.base $W/out.$w.$b > $W/diff.$w.$b && echo same || echo DIFFER)
    s=""; [ -f $W/out.$w.base/s.sam.zst ] && { cmp -s <(zstdcat $W/out.$w.base/s.sam.zst) <(zstdcat $W/out.$w.$b/s.sam.zst) && s="SAM text same" || s="SAM text DIFFER"; }
    echo "$w $b vs base: outputs $r ($(find $W/out.$w.base -type f | wc -l) files) $s"
  done
done
echo ONEBCG DONE
