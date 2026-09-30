#!/bin/bash
# T5b: kill -9 protal itself (not a timeout wrapper) during alignment; the rerun must realign.
set -u
source $(dirname "$0")/lib.sh
T=$W/t5; cd $T
I=$T/in
for fmt in zst gz sam; do
  d=kill2_$fmt; rm -rf $d; mkdir -p $d
  $P --db $DB -1 $I/big_R1.fq -2 $I/big_R2.fq --prefix big -o $d -t 1 --no_qcmsa --no_strains --sam_format $fmt > $d.log 2>&1 &
  pid=$!
  for i in $(seq 1 600); do grep -q 'Start parallel execution' $d.log 2>/dev/null && break; sleep 0.1; done
  sleep 2
  kill -9 $pid; wait $pid 2>/dev/null
  sleep 0.5; pgrep -f "prefix big -o $d" && echo "STILL RUNNING"
  echo "$fmt killed during alignment; files: $(cd $d && ls -la | awk 'NR>1{print $NF":"$5}' | grep -v '^\.' | tr '\n' ' ')"
  prun $d.re.log --db $DB -1 $I/big_R1.fq -2 $I/big_R2.fq --prefix big -o $d -t 2 --no_qcmsa --no_strains --sam_format $fmt
  grep -a -E 'Skip|All alignments|Align the|\[ERROR\]' $d.re.log | head -3
  s=$(ls $d/big.sam* | grep -v -e err -e partial)
  echo "after rerun: $(cd $d && ls | tr '\n' ' ') records=$(samtext $s | grep -vc '^@')"
done
