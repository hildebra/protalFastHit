#!/bin/bash
# Peak summed RSS (MB) of all protal processes while the command runs: memwatch.sh OUT command...
out=$1; shift
"$@" &
pid=$!
peak=0
while kill -0 $pid 2>/dev/null; do
  s=$(ps -C protal -o rss= | awk '{s+=$1} END {print int(s/1024)}')
  [ "${s:-0}" -gt "$peak" ] && peak=$s
  sleep 0.2
done
wait $pid; rc=$?
echo "peak protal RSS ${peak} MB, exit $rc" > $out
