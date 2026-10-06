#!/usr/bin/env bash
# summary.sh: one line per local run (set, threads, round): wall, user, max RSS, and the main stage timers.
cd $HOME/perf6
printf "%-4s %2s %2s %8s %8s %7s %7s %7s %7s %7s %7s\n" set t r wall user rss_gb index align prof seedfind alignh
for f in o.*_t*.r*.log; do
  b=${f%.log}; set=$(echo $b | sed -E 's/o\.([a-z]+)_t.*/\1/'); t=$(echo $b | sed -E 's/.*_t([0-9]+)\..*/\1/'); r=$(echo $b | sed -E 's/.*\.r([0-9]+)/\1/')
  wall=$(grep Elapsed $b.time | awk '{print $NF}'); user=$(grep 'User time' $b.time | awk '{print $NF}'); rss=$(grep 'Maximum resident' $b.time | awk '{printf "%.2f", $NF/1e6}')
  st() { bash $HOME/perf6/stages.sh $f | grep -m1 "^$1 " | awk '{print $NF}'; }
  printf "%-4s %2s %2s %8s %8s %7s %7s %7s %7s %7s %7s\n" $set $t $r $wall $user $rss "$(st 'Load Index')" "$(st 'Aligning reads')" "$(st 'Profiling')" "$(st 'Seed- and Anchor-finding')" "$(st 'Alignment handler')"
done
