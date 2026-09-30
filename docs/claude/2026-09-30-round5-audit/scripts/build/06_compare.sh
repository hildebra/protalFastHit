#!/bin/bash
# Are the baseline, AVX2 and static builds' outputs byte-identical on the mini DB?
# Three simulated samples (strains shared across samples, so strain MSAs and qcmsa run), each
# binary with -t 1 and -t 2; the .zst outputs are compared decompressed.
set -u
A=~/audit6/build
S=$A/src
C=$A/compare
rm -rf $C; mkdir -p $C/reads
cp -r ~/strain-build/mini_db/protal_db $C/db
G=~/strain-build/mini_db/gtdb_r226/simulation/genomes.tsv
cut -f1 $G | tail -n +2 > $C/accessions.txt; cat $C/accessions.txt | tr '\n' ' '; echo
mapfile -t acc < $C/accessions.txt
mk() { # name seed a1 w1 a2 w2 a3 w3
  printf 'accession\trelative_abundance\n%s\t%s\n%s\t%s\n%s\t%s\n' $3 $4 $5 $6 $7 $8 > $C/reads/$1.community.tsv
  python3 $S/scripts/mini_db/simulate_reads.py --genomes $G --community $C/reads/$1.community.tsv \
      --out_prefix $C/reads/$1 --pairs 20000 --seed $2 2>> $C/reads/sim.log
}
mk s1 7 ${acc[1]} 0.5 ${acc[3]} 0.3 ${acc[8]} 0.2
mk s2 8 ${acc[2]} 0.4 ${acc[3]} 0.3 ${acc[7]} 0.3
mk s3 9 ${acc[1]} 0.3 ${acc[4]} 0.4 ${acc[8]} 0.3
ls -la $C/reads | head
R1=$C/reads/s1_R1.fq,$C/reads/s2_R1.fq,$C/reads/s3_R1.fq
R2=$C/reads/s1_R2.fq,$C/reads/s2_R2.fq,$C/reads/s3_R2.fq
declare -A BIN=( [baseline]=$S/build/protal [avx2]=$S/build/protal_avx2 [static]=$S/build/protal_0.6.0_static [prebuilt]=$HOME/strain-build/bin/protal )
for b in baseline avx2 static prebuilt; do
  [ -x "${BIN[$b]}" ] || { echo "skip $b (no binary)"; continue; }
  for t in 1 2; do
    o=$C/out_${b}_t$t
    st=$(date +%s.%N)
    ( cd $C && taskset -c 0,1 "${BIN[$b]}" --db $C/db -1 $R1 -2 $R2 --prefix s1,s2,s3 -o $o -t $t \
        --qcmsa_script $S/scripts/qcmsa.py > $C/log_${b}_t$t.txt 2>&1 ); rc=$?
    echo "$b t$t rc=$rc $(echo "$(date +%s.%N) - $st" | bc 2>/dev/null) s; $(tail -1 $C/log_${b}_t$t.txt)"
    ( cd $o && find . -type f | sort | while read f; do
        case "$f" in *.zst) h=$(zstd -dc "$f" | md5sum | cut -d' ' -f1); n="${f%.zst} (zst, content)";; *) h=$(md5sum < "$f" | cut -d' ' -f1); n="$f";; esac
        echo "$h  $n"; done ) > $C/md5_${b}_t$t.txt
  done
done
echo "== files per run"; wc -l $C/md5_*.txt
cmpr() { if [ -f $C/md5_$1.txt ] && [ -f $C/md5_$2.txt ]; then
  d=$(diff $C/md5_$1.txt $C/md5_$2.txt | grep -c '^[<>]'); echo "$1 vs $2: $d differing lines"; diff $C/md5_$1.txt $C/md5_$2.txt | head -20; fi; }
cmpr baseline_t1 avx2_t1
cmpr baseline_t2 avx2_t2
cmpr baseline_t1 baseline_t2
cmpr avx2_t1 avx2_t2
cmpr baseline_t1 static_t1
cmpr baseline_t2 static_t2
cmpr baseline_t1 prebuilt_t1
echo "== profiles (baseline t1)"; cat $C/out_baseline_t1/s1.profile 2>/dev/null | head -8
echo DONE
