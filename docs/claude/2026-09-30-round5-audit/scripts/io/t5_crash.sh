#!/bin/bash
# T5: crash safety of SAM output (.partial + rename), rerun after a crash, write failures (EFBIG).
set -u
source $(dirname "$0")/lib.sh
T=$W/t5; rm -rf $T; mkdir -p $T/in; cd $T
I=$T/in
# a larger sample: sa repeated 40 times with unique names
for x in 1 2; do
  for k in $(seq 1 40); do awk -v k=$k 'NR%4==1{sub(/^@/, "@c" k "_")}{print}' $R/sa_R$x.fq; done > $I/big_R$x.fq
done
wc -l $I/big_R1.fq
for fmt in zst gz sam; do
  d=kill_$fmt; mkdir -p $d
  timeout 600 $P --db $DB -1 $I/big_R1.fq -2 $I/big_R2.fq --prefix big -o $d -t 1 --no_qcmsa --no_strains --sam_format $fmt > $d.log 2>&1 &
  pid=$!
  for i in $(seq 1 600); do grep -q 'Start parallel execution' $d.log 2>/dev/null && break; sleep 0.1; done
  sleep 1.5
  kill -9 $pid; wait $pid 2>/dev/null
  echo "$fmt killed during alignment; files: $(cd $d && ls -la | awk 'NR>1{print $NF":"$5}' | grep -v '^\.' | tr '\n' ' ')"
  prun $d.re.log --db $DB -1 $I/big_R1.fq -2 $I/big_R2.fq --prefix big -o $d -t 2 --no_qcmsa --no_strains --sam_format $fmt
  grep -a -E 'Skip|All alignments|Align the' $d.re.log | head -2
  s=$(ls $d/big.sam* | grep -v -e err -e partial)
  echo "after rerun: $(cd $d && ls | tr '\n' ' ') records=$(samtext $s | grep -vc '^@')"
done
echo "--- write failure: file size limit (EFBIG), each format, and with --full_sam_header"
for fmt in zst gz sam; do
  for extra in "" "--full_sam_header"; do
    d=efbig_${fmt}${extra:+_full}; mkdir -p $d
    bash -c "trap '' XFSZ; ulimit -f 200; exec $P --db $DB -1 $I/big_R1.fq -2 $I/big_R2.fq --prefix big -o $d -t 2 --no_qcmsa --no_strains --sam_format $fmt $extra" > $d.log 2>&1
    echo "$d rc=$? files: $(cd $d && ls | tr '\n' ' ') | $(grep -a -E '\[ERROR\]' $d.log | head -2 | tr '\n' '|' | cut -c1-250)"
  done
done
echo "--- unwritable output directory"
mkdir -p ro && chmod 555 ro
prun ro.log --db $DB -1 $R/sa_R1.fq -2 $R/sa_R2.fq --prefix sa -o ro -t 2 --no_qcmsa --no_strains
grep -a -E '\[ERROR\]|Cannot' ro.log | head -3
chmod 755 ro
