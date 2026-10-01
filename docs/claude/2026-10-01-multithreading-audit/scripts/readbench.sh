#!/usr/bin/env bash
# The input path's ceiling: readbench on the 5M-pair sample, gzip (BGZF, as the simulator writes it) and
# the same reads uncompressed, 1-6 reader threads, twice, niced. MTAUDIT lines per run in rb_<...>.err.
set -uo pipefail
M=${MT_AUDIT:-$HOME/mt-audit}; B=$HOME/bench071
R1=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz
R2=${R1%_R1.fq.gz}_R2.fq.gz
mkdir -p $M/plain $M/rb
[ -s $M/plain/R1.fq ] || zcat $R1 > $M/plain/R1.fq
[ -s $M/plain/R2.fq ] || zcat $R2 > $M/plain/R2.fq
head -c 4 $R1 | od -An -tx1
echo "bgzf extra field: $(head -c 16 $R1 | od -An -tx1 | tr -s ' ')"
cat $M/plain/R1.fq $M/plain/R2.fq > /dev/null; cat $R1 $R2 > /dev/null   # into the page cache
for rep in 1 2; do
  for kind in gz plain; do
    for t in 1 2 3 4 6; do
      if [ $kind = gz ]; then a=$R1; b=$R2; else a=$M/plain/R1.fq; b=$M/plain/R2.fq; fi
      out=$(nice -n 5 $M/readbench/readbench $a $b $t 2> $M/rb/rb_${kind}_${t}_$rep.err)
      echo -e "$kind\trep$rep\t$out\twait_s\t$(awk -F'\t' '$2=="reader_wait"{print $3}' $M/rb/rb_${kind}_${t}_$rep.err)\thold_s\t$(awk -F'\t' '$2=="reader_hold"{print $3}' $M/rb/rb_${kind}_${t}_$rep.err)\tinflate_busy_s\t$(awk -F'\t' '$2=="inflate_busy"{print $3}' $M/rb/rb_${kind}_${t}_$rep.err)\tunderflow_s\t$(awk -F'\t' '$2=="underflow_wait"{print $3}' $M/rb/rb_${kind}_${t}_$rep.err)\tload\t$(cut -d' ' -f1 /proc/loadavg)"
    done
  done
done
