#!/bin/bash
# A poor man's wall-clock sampler (no perf here): run the profiling build (debug info) under gdb and
# interrupt it N times, recording the innermost frames, inlined ones included. One thread, pinned to one vCPU.
#   sample.sh LABEL DB READS_DIR SAMPLES [DELAY_S] [INTERVAL_S]
# Sampling starts DELAY_S (default 6, after the index load) into the run. Output: $PERF_DIR/sample/LABEL/{raw.txt,top.txt}.
source "$(dirname "$0")/env.sh"
label=$1; db=$2; reads=$3; n=$4; delay=${5:-6}; interval=${6:-0.05}
out=$PERF_DIR/sample/$label; rm -rf $out; mkdir -p $out/o
r1=$(ls $reads/*_R1.fq* | head -1); r2=$(ls $reads/*_R2.fq* | head -1)
cat > $out/gdb.cmd <<GDB
set pagination off
set confirm off
handle SIGINT stop nopass
run
set \$i = 0
while \$i < $n
  echo ===\n
  thread 1
  bt 8
  set \$i = \$i + 1
  continue
end
kill
quit
GDB
( sleep $delay; for i in $(seq 1 $n); do pid=$(pgrep -x protal_avx2 | head -1); [ -n "$pid" ] && kill -INT $pid; sleep $interval; done ) &
gdb -batch -x $out/gdb.cmd --args taskset -c 2 $PROF_BIN --db $db -1 $r1 -2 $r2 -o $out/o/out -t 1 --no_qcmsa --no_profile > $out/raw.txt 2>/dev/null
wait
bash $here/sample_report.sh $out/raw.txt $out

