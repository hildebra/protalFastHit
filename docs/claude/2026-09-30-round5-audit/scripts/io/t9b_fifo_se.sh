#!/bin/bash
# T9b: single-end FIFO input, and a process substitution as protal users write it.
set -u
source $(dirname "$0")/lib.sh
T=$W/t9; cd $T; I=$T/in
rm -f $I/s1; mkfifo $I/s1
( cat $R/sa_R1.fq > $I/s1 ) &
prun fifose.log --db $DB -1 $I/s1 --model_se $DB/model_pe.xml --prefix fifose -o out -t 2 --no_qcmsa --no_strains --sam_format sam
kill %1 2>/dev/null; wait 2>/dev/null
grep -a -E 'unrecog|malformed|\[ERROR\]' fifose.log | head -3
gzip -c $R/sa_R1.fq > $I/a1.gz; gzip -c $R/sa_R2.fq > $I/a2.gz
timeout 600 $P --db $DB -1 <(zcat $I/a1.gz) -2 <(zcat $I/a2.gz) --prefix procsub -o out -t 2 --no_qcmsa --no_strains --sam_format sam > procsub.log 2>&1
echo "process substitution rc=$? records=$(grep -vc '^@' out/procsub.sam 2>/dev/null)"
grep -a -E 'unrecog|malformed|\[ERROR\]' procsub.log | head -3
