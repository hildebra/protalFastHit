#!/bin/bash
# E5: malformed internal_taxonomy.dmp variants of the unpacked mini DB, protal on 4000 pairs
# (every 10th pair of E1 sample s_1, so that all its genomes contribute).
# Run: nohup bash e5_tax.sh > ~/audit6/simulator/e5/e5.out 2>&1 &
set -u
W=~/audit6/simulator
cd $W/e5
zcat $W/e1/sim/reads/s_1_R1.fq.gz | awk '{r=int((NR-1)/4)} r%10==0' > r1.fq
zcat $W/e1/sim/reads/s_1_R2.fq.gz | awk '{r=int((NR-1)/4)} r%10==0' > r2.fq
echo "pairs: $(( $(wc -l < r1.fq) / 4 ))"
B=$W/e5/dbsrc
mk() {  # mk NAME SED-SCRIPT
  rm -rf db_$1; mkdir db_$1
  for f in index.prx.zst model_pe.xml reference.fna reference.map unique_kmers.tsv; do ln -s $B/$f db_$1/$f; done
  sed "$2" $B/internal_taxonomy.dmp > db_$1/internal_taxonomy.dmp
}
mk base ''
mk cycle 's/^13\t11\t/13\t15\t/'
mk noroot 's/^4\t4\t/4\t2\t/'
mk rootgone '/^4\t4\t/d'
mk slash 's/s__Mockella alpha/s__Mockella al\/pha/'
mk semicolon 's/s__Mockella alpha/s__Mockella al;pha/'
mk dupname 's/s__Mockella beta/s__Mockella alpha/'
mk rankcase 's/^15\t13\t0\tg__Mockella\tgenus/15\t13\t0\tg__Mockella\tGenus/'
for v in ${VARIANTS:-base cycle noroot rootgone slash semicolon dupname rankcase}; do
  echo "=================== $v"
  diff $B/internal_taxonomy.dmp db_$v/internal_taxonomy.dmp | grep '^[<>]'
  rm -rf out_$v
  timeout -s KILL 90 ~/strain-build/bin/protal --db db_$v -1 r1.fq -2 r2.fq -o out_$v -t 2 --no_strains --profile_truth truth.txt > log_$v.txt 2> err_$v.txt
  echo "exit $?"
  grep -h "Invalid taxonomy\|TP \|Species unknown\|failed\|rror" log_$v.txt err_$v.txt | sort -u | head -8
  echo "-- profile:"; cat out_$v/*.profile 2>/dev/null
  echo "-- misc statistics files:"; (cd out_$v/misc 2>/dev/null && find . -name "*.statistics.tsv" | sort)
  echo "-- truth_annotated (truth prediction taxon name):"; cut -f1,2,4,5 out_$v/*.truth_annotated 2>/dev/null | tail -n +2
done
echo DONE
