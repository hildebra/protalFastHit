#!/bin/bash
# Thread scaling of the alignment stage: threads.sh BIN DB READS_DIR REPS "1 2 3 4 5 6"
# Median wall, user CPU and loop time per thread count; a fixed part (wall - loop) shows the serial start-up.
source "$(dirname "$0")/env.sh"
bin=$1; db=$2; reads=$3; reps=$4; ts=${5:-"1 2 3 4 5 6"}
r1=$(ls $reads/*_R1.fq* | head -1); r2=$(ls $reads/*_R2.fq* | head -1)
d=$PERF_DIR/runs2/threads; mkdir -p $d
for rep in $(seq 1 $reps); do
  for t in $ts; do
    rm -rf $d/out*
    /usr/bin/time -f "$t %e %U %S" -o $d/time.txt $bin --db $db -1 $r1 -2 $r2 -o $d/out -t $t --no_qcmsa --no_profile --verbose > $d/stdout.log 2> $d/stderr.log
    loop=$(grep "^OMP Loop handler took" $d/stdout.log | sed -E 's/.* took //; s/ms//; s/s / /' | awk '{ if (NF == 2) print $1 + $2 / 1000; else print $1 / 1000 }')
    echo "$(cat $d/time.txt) $loop"
  done
done > $d/all.txt
echo "threads  wall  user  sys  loop   (medians of $reps)"
for t in $ts; do
  awk -v t=$t '$1 == t { w[++n] = $2; u[n] = $3; s[n] = $4; l[n] = $5 } END {
    asort(w); asort(u); asort(s); asort(l); m = int((n + 1) / 2); printf "%5d %7.2f %6.2f %5.2f %6.2f\n", t, w[m], u[m], s[m], l[m] }' $d/all.txt
done
rm -rf $d/out*
