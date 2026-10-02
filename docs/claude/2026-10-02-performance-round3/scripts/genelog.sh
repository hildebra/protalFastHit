#!/usr/bin/env bash
a=$HOME/mt-work/perf3/pe5.ref/s.profile.gene.log; b=$HOME/mt-work/perf3/pe5.lr/s.profile.gene.log
diff $a $b | head -8; wc -l $a
echo "sorted diff lines: $(diff <(sort $a) <(sort $b) | wc -l)"
head -3 $a | cut -c1-200
# the reference run against itself: run the reference again at -t 6 and compare
W=$HOME/mt-work/perf3; P=$HOME/bench071/samples/points; R5=$P/rl150_p500000/sim/reads/rl150_p500000_s_1
rm -rf $W/pe5.ref2; nice $W/ref/build/protal --db $HOME/bench071/V071/protal_db -1 ${R5}_R1.fq.gz -2 ${R5}_R2.fq.gz --prefix s -o $W/pe5.ref2 -t 6 --no_qcmsa > /dev/null 2>&1
cmp -s $a $W/pe5.ref2/s.profile.gene.log && echo "reference twice: gene log same" || echo "reference twice: gene log DIFFERS ($(diff <(sort $a) <(sort $W/pe5.ref2/s.profile.gene.log) | wc -l) sorted diff lines)"
