#!/usr/bin/env bash
# readbench against the work tree (parallel BGZF): BGZF input at 1, 2 and 3 inflating threads per file
# (PROTAL_INFLATE_THREADS), uncompressed and one-member gzip for reference; 2 consumers, then 3 with no work; 4 cores.
G=$HOME/gzpar; D=$G/data; T="taskset -c 0-3 nice -n 5"
S=$G/work
cd $G/rb
$T g++ -std=c++20 -O3 -march=x86-64-v2 -fopenmp -I$S/src -I$S/src/IO -I$S/src/SequenceUtils -I$S/src/Utilities -I$S/lib \
   readbench.cpp $S/src/IO/FastxReader.cpp -lisal -lzstd -lpthread -o readbench_work 2> build_work.err || { head -30 build_work.err; exit 1; }
cat $D/R1.fq $D/R2.fq $D/R1.bgzf.gz $D/R2.bgzf.gz $D/R1.member.gz $D/R2.member.gz > /dev/null
uptime
for rep in 1 2 3; do
  for t in 2 3; do
    for case in fq member1 bgzf1 bgzf2 bgzf3; do
      case $case in
        fq) k=fq; n=1 ;; member1) k=member.gz; n=1 ;; bgzf1) k=bgzf.gz; n=1 ;; bgzf2) k=bgzf.gz; n=2 ;; bgzf3) k=bgzf.gz; n=3 ;;
      esac
      out=$(PROTAL_INFLATE_THREADS=$n $T ./readbench_work $D/R1.$k $D/R2.$k $t 2>/dev/null)
      echo -e "$case\tconsumers $t\trep$rep\t$(echo "$out" | awk -F'\t' '{print "pairs/s", $10, "s", $8, "pairs", $6}')"
    done
  done
done
uptime
