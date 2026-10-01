#!/bin/bash
# The experiment builds, one after the other (each takes 4-5 minutes, most of it the serial LTO link):
# flank counting, the cheap fixes, then PGO. Progress goes to $PERF_DIR/experiments.out.
source "$(dirname "$0")/env.sh"
log=$PERF_DIR/experiments.out
for b in build_flankcount build_quick build_pgo; do
  echo "$(date +%T) $b start" >> $log
  bash $here/$b.sh > $PERF_DIR/$b.out 2>&1
  echo "$(date +%T) $b finished rc=$?" >> $log
done
echo "ALL DONE" >> $log
