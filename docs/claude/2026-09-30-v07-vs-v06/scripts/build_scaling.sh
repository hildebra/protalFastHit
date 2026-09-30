#!/usr/bin/env bash
# More database builds of the tuning world (after build.sh tune): each version on 1 thread, and 0.7's
# single file at zstd level 3 (the training database's level) on all threads. Timed as build.sh's.
set -euo pipefail
B=${BENCH:-$HOME/bench07}
P7=${P7:-$HOME/fix-build/bin/protal}
P6=${P6:-$HOME/protal-0.6.0a/src/build/protal}
T=${T:-6}
W=$B/tune
tm() { local step=$1; shift; echo "$(date +%T) $step"; /usr/bin/time -v -o $B/logs/tune.$step.time "$@" > $B/logs/tune.$step.log 2>&1; }
one() {  # step, binary, threads, extra options
  local step=$1 p=$2 t=$3; shift 3
  rm -rf $W/$step; cp -r $W/conv $W/$step
  [ $p = $P6 ] && mv $W/$step/model_pe.xml $W/$step/model.xml
  tm $step $p --build --no_profile -t $t --db $W/$step --reference $W/$step/reference.fna --full_reference $W/$step/full_reference.fna "$@"
  du -sb $W/$step >> $B/logs/tune.sizes.txt
  rm -rf $W/$step
}
one build07raw_t1 $P7 1 --no_compress
one build06_t1 $P6 1
one build07_l3 $P7 $T --compress_level 3
