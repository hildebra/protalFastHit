#!/usr/bin/env bash
# Today's reader (HEAD): pairs per second from the input path alone, 2 consumer threads (+ 2 inflating), 4 cores.
G=$HOME/gzpar; D=$G/data; RB=${RB:-$G/rb/readbench}; T="taskset -c 0-3 nice -n 5"
cat $D/R1.fq $D/R2.fq $D/R1.bgzf.gz $D/R2.bgzf.gz $D/R1.member.gz $D/R2.member.gz > /dev/null
uptime
for rep in 1 2; do
  for kind in fq bgzf.gz member.gz; do
    for t in 1 2; do
      echo -e "$kind\trep$rep\t$($T $RB $D/R1.$kind $D/R2.$kind $t 2>/dev/null)"
    done
  done
done
uptime
