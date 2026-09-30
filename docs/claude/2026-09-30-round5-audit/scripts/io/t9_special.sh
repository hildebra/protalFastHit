#!/bin/bash
# T9: special input files: named pipes (FIFOs), unreadable files, a directory; a 3 MB read past the
# long-read check (reads longer than the 1 MB inflate block, plain/gzip/BGZF).
set -u
source $(dirname "$0")/lib.sh
T=$W/t9; rm -rf $T; mkdir -p $T/in; cd $T
I=$T/in; A1=$R/sa_R1.fq; A2=$R/sa_R2.fq
BASE=$W/t1/sam_t1.rec
res() { # res NAME
  local name=$1 s=out/$1.sam rec=none same=-
  [ -e $s ] && { records $s > $name.rec; rec=$(wc -l < $name.rec); same=$(cmp -s $name.rec $BASE && echo SAME || echo diff); }
  printf '%-10s %s records=%s vs_base=%s | %s\n' $name "$(cat $name.rc)" $rec $same "$(grep -a -E '\[ERROR\]|rror|malformed|unrecog|differ' $name.log | head -3 | tr '\n' '|' | cut -c1-300)"
}
echo "--- FIFOs (plain and gzip data)"
mkfifo $I/p1 $I/p2
( cat $A1 > $I/p1 ) & ( cat $A2 > $I/p2 ) &
prun fifo.log --db $DB -1 $I/p1 -2 $I/p2 --prefix fifo -o out -t 2 --no_qcmsa --no_strains --sam_format sam > fifo.rc
kill %1 %2 2>/dev/null; wait 2>/dev/null
res fifo
rm -f $I/p1 $I/p2; mkfifo $I/p1 $I/p2
( gzip -c $A1 > $I/p1 ) & ( gzip -c $A2 > $I/p2 ) &
prun fifogz.log --db $DB -1 $I/p1 -2 $I/p2 --prefix fifogz -o out -t 2 --no_qcmsa --no_strains --sam_format sam > fifogz.rc
kill %1 %2 2>/dev/null; wait 2>/dev/null
res fifogz
echo "--- unreadable input files"
cp $A1 $I/nr_R1.fq; cp $A2 $I/nr_R2.fq; chmod 000 $I/nr_R1.fq $I/nr_R2.fq
prun nr.log --db $DB -1 $I/nr_R1.fq -2 $I/nr_R2.fq --prefix nr -o out -t 2 --no_qcmsa --no_strains --sam_format sam > nr.rc; res nr
prun nrse.log --db $DB -1 $I/nr_R1.fq --model_se $DB/model_pe.xml --prefix nrse -o out -t 2 --no_qcmsa --no_strains --sam_format sam > nrse.rc; res nrse
gzip -c $A1 > $I/nrgz_R1.fq.gz; gzip -c $A2 > $I/nrgz_R2.fq.gz; chmod 000 $I/nrgz_R1.fq.gz $I/nrgz_R2.fq.gz
prun nrgz.log --db $DB -1 $I/nrgz_R1.fq.gz -2 $I/nrgz_R2.fq.gz --prefix nrgz -o out -t 2 --no_qcmsa --no_strains --sam_format sam > nrgz.rc; res nrgz
python3 $S/mkreads.py bgzf $A1 $I/nrb_R1.fq.gz; python3 $S/mkreads.py bgzf $A2 $I/nrb_R2.fq.gz; chmod 000 $I/nrb_R1.fq.gz $I/nrb_R2.fq.gz
prun nrb.log --db $DB -1 $I/nrb_R1.fq.gz -2 $I/nrb_R2.fq.gz --prefix nrb -o out -t 2 --no_qcmsa --no_strains --sam_format sam > nrb.rc; res nrb
chmod 644 $I/nr*
echo "--- a directory as input"
mkdir -p $I/d1 $I/d2
prun dir.log --db $DB -1 $I/d1 -2 $I/d2 --prefix dir -o out -t 2 --no_qcmsa --no_strains --sam_format sam > dir.rc; res dir
echo "--- a 3 MB read after the first 100 reads (the long-read check looks at 100)"
python3 - <<'EOF' > $I/long.fq
import random
r = random.Random(5)
src = open("/home/falk/audit6/io/reads/sa_R1.fq").read().split("\n")
recs = ["\n".join(src[i:i+4]) for i in range(0, 1200, 4)]
L = 3 * 1024 * 1024
s = "".join(r.choice("ACGT") for _ in range(L))
out = recs[:150] + ["@long\n" + s + "\n+\n" + "I" * L] + recs[150:]
print("\n".join(out))
EOF
head -n 1200 $A1 > $I/short.fq
gzip -c $I/long.fq > $I/long.fq.gz; python3 $S/mkreads.py bgzf $I/long.fq $I/longb.fq.gz
for v in short long; do
  prun $v.log --db $DB -1 $I/$v.fq --model_se $DB/model_pe.xml --prefix $v -o out -t 2 --no_qcmsa --no_strains --sam_format sam > $v.rc
done
for v in longgz:long.fq.gz longb:longb.fq.gz; do
  n=${v%%:*}; f=${v#*:}
  prun $n.log --db $DB -1 $I/$f --model_se $DB/model_pe.xml --prefix $n -o out -t 2 --no_qcmsa --no_strains --sam_format sam > $n.rc
done
for v in short long longgz longb; do
  s=out/$v.sam
  echo "$v $(cat $v.rc) qnames=$( [ -e $s ] && grep -v '^@' $s | cut -f1 | sort -u | wc -l) long_records=$( [ -e $s ] && grep -c '^long' $s) shortreads_same=$( [ -e $s ] && diff <(grep -v '^@' $s | grep -v '^long' | LC_ALL=C sort) <(grep -v '^@' out/short.sam | LC_ALL=C sort) >/dev/null && echo yes || echo no) | $(grep -a -E '\[ERROR\]|rror' $v.log | head -2 | tr '\n' '|')"
done
