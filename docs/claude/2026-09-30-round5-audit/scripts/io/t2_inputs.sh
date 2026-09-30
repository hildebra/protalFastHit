#!/bin/bash
# T2: FASTQ/FASTA input variants: line endings, compression, truncation, corruption, pairing.
set -u
source $(dirname "$0")/lib.sh
T=$W/t2; rm -rf $T; mkdir -p $T/in; cd $T
I=$T/in
A1=$R/sa_R1.fq; A2=$R/sa_R2.fq
BASE=$W/t1/sam_t1.rec   # records of the plain run, sorted
n=$(wc -l < $A1); half=$(( (n/8)*4 ))   # a whole number of records

mk() { # mk NAME R1 R2 : both files prepared by caller
  :
}
# ---- variants (R1 and R2 transformed alike unless noted)
for x in 1 2; do
  src=$R/sa_R$x.fq
  sed 's/$/\r/' $src > $I/crlf_R$x.fq
  head -c -1 $src > $I/nonl_R$x.fq                       # no trailing newline
  gzip -c $src > $I/gz_R$x.fq.gz
  { head -n $half $src | gzip -c; tail -n +$((half+1)) $src | gzip -c; } > $I/multi_R$x.fq.gz
  python3 $S/mkreads.py bgzf $src $I/bgzf_R$x.fq.gz
  { head -n $half $src > $I/h$x; tail -n +$((half+1)) $src > $I/t$x; }
  python3 $S/mkreads.py bgzf $I/h$x $I/hb$x; gzip -c $I/t$x > $I/tg$x
  cat $I/hb$x $I/tg$x > $I/bgzfplusgz_R$x.fq.gz             # BGZF members, then a plain gzip member
  cat $I/tg$x $I/hb$x > $I/gzplusbgzf_R$x.fq.gz             # plain gzip first, then BGZF
  python3 $S/mkreads.py bgzf $I/h$x $I/hb2$x; python3 $S/mkreads.py bgzf $I/t$x $I/tb2$x
  cat $I/hb2$x $I/tb2$x > $I/bgzf2_R$x.fq.gz                # two BGZF files (EOF block in the middle)
  awk 'NR%4==2{print tolower($0);next}{print}' $src > $I/lower_R$x.fq
  awk 'BEGIN{srand(7)} NR%4==2{s=$0; for(i=5;i<=length(s);i+=17){s=substr(s,1,i-1) "N" substr(s,i+1)}; print s; next}{print}' $src > $I/nbase_R$x.fq
  awk 'NR%4==1{print ">" substr($0,2)} NR%4==2{print}' $src > $I/fasta_R$x.fa
  awk 'NR%4==1{print ">" substr($0,2)} NR%4==2{print substr($0,1,50); print substr($0,51)}' $src > $I/fastaml_R$x.fa  # multi-line FASTA
  cp $I/gz_R$x.fq.gz $I/trunc_R$x.fq.gz
done
# truncations and corruptions (R1 only is damaged; R2 fine)
sz=$(stat -c %s $I/gz_R1.fq.gz); head -c $((sz/2)) $I/gz_R1.fq.gz > $I/trunc_R1.fq.gz
head -c $((sz-4)) $I/gz_R1.fq.gz > $I/truncend_R1.fq.gz; cp $I/gz_R2.fq.gz $I/truncend_R2.fq.gz   # cut into the trailer
cp $I/gz_R1.fq.gz $I/corrupt_R1.fq.gz; printf '\x55' | dd of=$I/corrupt_R1.fq.gz bs=1 seek=$((sz/2)) conv=notrunc 2>/dev/null; cp $I/gz_R2.fq.gz $I/corrupt_R2.fq.gz
cp $I/gz_R1.fq.gz $I/crc_R1.fq.gz; printf '\x00\x00\x00\x00' | dd of=$I/crc_R1.fq.gz bs=1 seek=$((sz-8)) conv=notrunc 2>/dev/null; cp $I/gz_R2.fq.gz $I/crc_R2.fq.gz  # CRC zeroed
bsz=$(stat -c %s $I/bgzf_R1.fq.gz)
head -c $((bsz/2)) $I/bgzf_R1.fq.gz > $I/btrunc_R1.fq.gz; cp $I/bgzf_R2.fq.gz $I/btrunc_R2.fq.gz
head -c $((bsz-28)) $I/bgzf_R1.fq.gz > $I/bnoeof_R1.fq.gz; cp $I/bgzf_R2.fq.gz $I/bnoeof_R2.fq.gz    # cut at a block boundary: EOF block missing
cp $I/bgzf_R1.fq.gz $I/bcrc_R1.fq.gz; printf '\x00\x00\x00\x00' | dd of=$I/bcrc_R1.fq.gz bs=1 seek=$((bsz-36)) conv=notrunc 2>/dev/null; cp $I/bgzf_R2.fq.gz $I/bcrc_R2.fq.gz
# multi-member gzip cut exactly at the member boundary: indistinguishable from a complete file
m1=$(head -n $half $A1 | gzip -c | wc -c); head -c $m1 $I/multi_R1.fq.gz > $I/mcut_R1.fq.gz; cp $I/multi_R2.fq.gz $I/mcut_R2.fq.gz
{ cat $I/gz_R1.fq.gz; printf 'garbage after the gzip member\n'; } > $I/trail_R1.fq.gz; cp $I/gz_R2.fq.gz $I/trail_R2.fq.gz
# pairing problems
cp $A1 $I/r2short_R1.fq; head -n $((n-4)) $A2 > $I/r2short_R2.fq            # R2 one read short
cp $A1 $I/r2short33_R1.fq; head -n $((n-132)) $A2 > $I/r2short33_R2.fq      # R2 33 reads short
cp $A1 $I/r1long_R1.fq; head -n $((n-128)) $A2 > $I/r1long_R2.fq            # 32 short: a whole batch
awk 'NR%4==1{sub(/\/2$/,""); print "@x" substr($0,2) "/2"; next}{print}' $A2 > $I/names_R2.fq; cp $A1 $I/names_R1.fq   # R2 names differ
# shuffled R2 (pairs broken, all names mismatched)
paste - - - - < $A2 | awk 'BEGIN{srand(3)}{print rand() "\t" $0}' | sort -k1,1 | cut -f2- | tr '\t' '\n' > $I/shuf_R2.fq; cp $A1 $I/shuf_R1.fq
# empty inputs
: > $I/empty_R1.fq; : > $I/empty_R2.fq
: > $I/empty1_R1.fq; cp $A2 $I/empty1_R2.fq
printf '' | gzip -c > $I/emptygz_R1.fq.gz; cp $I/emptygz_R1.fq.gz $I/emptygz_R2.fq.gz
# blank line in the middle (after record 10) / at the end
awk 'NR==41{print ""}{print}' $A1 > $I/blank_R1.fq; awk 'NR==41{print ""}{print}' $A2 > $I/blank_R2.fq
{ cat $A1; echo; } > $I/blankend_R1.fq; { cat $A2; echo; } > $I/blankend_R2.fq
# a record whose quality line starts with '@' and a header "@" alone
# not a FASTQ at all
printf 'hello\nworld\n' > $I/notfq_R1.fq; cp $I/notfq_R1.fq $I/notfq_R2.fq

run() { # run NAME R1 R2 [extra]
  local name=$1 r1=$2 r2=$3; shift 3
  prun $name.log --db $DB -1 $r1 -2 $r2 --prefix $name -o out -t 2 --no_qcmsa --no_strains --sam_format sam "$@" > $name.rc
  local sam=out/$name.sam rec=0 same=-
  if [ -e $sam ]; then records $sam > $name.rec; rec=$(wc -l < $name.rec); same=$(cmp -s $name.rec $BASE && echo SAME || echo diff); fi
  printf '%-14s %s sam=%s records=%s vs_base=%s | %s\n' $name "$(cat $name.rc)" "$([ -e $sam ] && echo yes || echo no)" $rec $same \
     "$(grep -i -E 'error|abort|truncat|corrupt|fail|differ|not the same|malformed|unrecognized|warn' $name.log | grep -v -E '^(snp|msa)' | head -3 | tr '\n' '|' | cut -c1-300)"
}
for v in crlf nonl; do run $v $I/${v}_R1.fq $I/${v}_R2.fq; done
for v in gz multi bgzf bgzfplusgz gzplusbgzf bgzf2 trunc truncend corrupt crc btrunc bnoeof bcrc mcut trail emptygz; do run $v $I/${v}_R1.fq.gz $I/${v}_R2.fq.gz; done
for v in lower nbase r2short r2short33 r1long names shuf empty empty1 blank blankend notfq; do run $v $I/${v}_R1.fq $I/${v}_R2.fq; done
for v in fasta fastaml; do run $v $I/${v}_R1.fa $I/${v}_R2.fa; done
echo "--- lowercase: first differing records"
diff <(cut -f1-9 lower.rec) <(cut -f1-9 $BASE) | head -4
diff lower.rec $BASE | head -4 | cut -c1-200
echo "--- fasta vs base (first diffs, cols 1-9)"
diff <(cut -f1-9 fasta.rec) <(cut -f1-9 $BASE) | head -4
