#!/bin/bash
# Re-time protal on the floor experiment's 150 Mb samples (12 ONT samples, 1.8 GB of reads) at the three ONT floors,
# three rounds in rotating order; outputs to ~/ontu/timing (the trace's own outputs stay). Cores 0-3, niced.
R=~/ontu/rw_scratch
T=~/ontu/trace/d150000000
for r in 1 2 3; do
  case $r in 1) order="0.85 0.80 0.75";; 2) order="0.80 0.75 0.85";; 3) order="0.75 0.85 0.80";; esac
  for A in $order; do
    O=~/ontu/timing/a$A
    rm -rf $O && mkdir -p $O
    sed "s#$T/out_a$A#$O#g" $T/out_a$A/samples.map > $O/samples.map
    /usr/bin/time -f "round $r -a $A: %e s wall, %U s user, %S s sys, %M KB" taskset -c 0-3 nice -n 5 ~/ontu/build/protal \
      --db $R/training_db --map $O/samples.map -t 4 --no_strains --no_qcmsa --sparse_index emulate -a $A > $O/protal.log 2>$O/time.txt
    tail -1 $O/time.txt; uptime | sed 's/.*load/  load/'
  done
done
