#!/bin/bash
# Start-up cost with and without the cPMML patches (patches/): 1000 pairs, 1 thread, alternated,
# plus instruction counts. The profiles must be identical.
#   model_test.sh DB READS_DIR
# Builds $PERF_DIR/build-exp from $PERF_DIR/src with the patches applied.
set -e
source "$(dirname "$0")/env.sh"
db=$1; reads=$2
rm -rf $PERF_DIR/src-exp; cp -a $PERF_DIR/src $PERF_DIR/src-exp
(cd $PERF_DIR/src-exp && for p in $here/../patches/*.patch; do patch -p1 < $p; done)
cmake -S $PERF_DIR/src-exp -B $PERF_DIR/build-exp -G Ninja -DCMAKE_BUILD_TYPE=Release > $PERF_DIR/build-exp.cmake.log 2>&1
cmake --build $PERF_DIR/build-exp --target protal_avx2 -j "$(nproc)" > $PERF_DIR/build-exp.log 2>&1
tiny=$PERF_DIR/reads/tiny; mkdir -p $tiny
zcat $reads/*_R1.fq.gz | head -4000 > $tiny/tiny_R1.fq
zcat $reads/*_R2.fq.gz | head -4000 > $tiny/tiny_R2.fq
set +e
for rep in 1 2 3; do
  for v in rel exp; do
    out=$PERF_DIR/runs/model_${v}_$rep; rm -rf $out; mkdir -p $out
    /usr/bin/time -f "$v rep$rep: %e s wall %U s user" $PERF_DIR/build-$v/protal_avx2 --db $db -1 $tiny/tiny_R1.fq -2 $tiny/tiny_R2.fq \
      -o $out/out -t 1 --no_qcmsa > $out/stdout.log 2> $out/time.txt
    tail -1 $out/time.txt
  done
done
cmp $PERF_DIR/runs/model_rel_1/out/tiny_R.profile $PERF_DIR/runs/model_exp_1/out/tiny_R.profile && echo "profiles identical"
for v in rel exp; do
  valgrind --tool=callgrind --callgrind-out-file=$PERF_DIR/runs/model_${v}_1/cg.out $PERF_DIR/build-$v/protal_avx2 --db $db \
    -1 $tiny/tiny_R1.fq -2 $tiny/tiny_R2.fq -o $PERF_DIR/runs/model_${v}_1/cgout -t 1 --no_qcmsa > /dev/null 2>&1
  echo "$v: $(callgrind_annotate --inclusive=yes --threshold=100 $PERF_DIR/runs/model_${v}_1/cg.out 2>/dev/null | grep -E 'PROGRAM TOTALS|__cxa_throw \[' | awk '{print $1}' | tr '\n' ' ')(total, __cxa_throw instructions)"
done
