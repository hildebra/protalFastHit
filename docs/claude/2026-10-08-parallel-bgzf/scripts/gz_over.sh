#!/usr/bin/env bash
# Consumers that also work (work_us per pair) on all 4 cores, as aligners do: does the input keep up when its inflating
# threads share the cores with busy consumers? 4 consumers + 2 inflating threads on 4 cores.
G=$HOME/gzpar; D=$G/data; RB=${RB:-$G/rb/readbench}; T="taskset -c 0-3 nice -n 5"
uptime
for rep in 1 2; do
  for w in 1.2 1.8 3.0; do
    for kind in fq member.gz bgzf.gz; do
      echo -e "$kind\trep$rep\t$($T $RB $D/R1.$kind $D/R2.$kind 4 $w 2>/dev/null)"
    done
  done
done
uptime
