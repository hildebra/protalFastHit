#!/usr/bin/env bash
# Standard gzip input (one member, zlib-ng path) against BGZF (libdeflate path): readbench with 1 and 2
# readers, 3 alternated reps, niced. Also how fast each inflates alone (inflate_busy per file).
set -uo pipefail
M=${MT_AUDIT:-$HOME/mt-audit}; B=$HOME/bench071
R1=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz
R2=${R1%_R1.fq.gz}_R2.fq.gz
cd $M/plain
[ -s R1.std.gz ] || nice pigz -p 2 -6 -c R1.fq > R1.std.gz
[ -s R2.std.gz ] || nice pigz -p 2 -6 -c R2.fq > R2.std.gz
ls -la R1.std.gz $R1
head -c 16 R1.std.gz | od -An -tx1
cat R1.std.gz R2.std.gz $R1 $R2 > /dev/null
for rep in 1 2 3; do
  for kind in bgzf std; do
    for t in 1 2; do
      if [ $kind = bgzf ]; then a=$R1; b=$R2; else a=$M/plain/R1.std.gz; b=$M/plain/R2.std.gz; fi
      e=$M/rb/rbg_${kind}_${t}_$rep.err
      out=$(nice -n 5 $M/readbench/readbench $a $b $t 2> $e)
      g() { awk -F'\t' -v k=$1 '$2==k{print $3}' $e; }
      echo -e "$kind\trep$rep\t$out\twait_s\t$(g reader_wait)\thold_s\t$(g reader_hold)\tinflate_busy_s\t$(g inflate_busy)\tunderflow_s\t$(g underflow_wait)\tload\t$(cut -d' ' -f1 /proc/loadavg)"
    done
  done
done
