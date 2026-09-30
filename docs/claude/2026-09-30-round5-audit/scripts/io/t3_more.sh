#!/bin/bash
# T3: single-end variants, damaged multi-member gzip, other compressors, long lines, lowercase,
# preload off, full header, zero-alignment SAMs.
set -u
source $(dirname "$0")/lib.sh
T=$W/t3; rm -rf $T; mkdir -p $T/in; cd $T
I=$T/in; I2=$W/t2/in
A1=$R/sa_R1.fq; A2=$R/sa_R2.fq
BASE=$W/t1/sam_t1.rec
n=$(wc -l < $A1); half=$(( (n/8)*4 ))
which bzip2 xz || true
# multi-member gzip whose 2nd member header is damaged (first byte 0x1f -> 0x00)
for x in 1 2; do
  m=$(head -n $half $R/sa_R$x.fq | gzip -c | wc -c)
  cp $I2/multi_R$x.fq.gz $I/mbad_R$x.fq.gz
  printf '\x00' | dd of=$I/mbad_R$x.fq.gz bs=1 seek=$m conv=notrunc 2>/dev/null
done
command -v bzip2 >/dev/null && { bzip2 -c $A1 > $I/bz_R1.fq.bz2; bzip2 -c $A2 > $I/bz_R2.fq.bz2; }
command -v xz >/dev/null && { xz -c $A1 > $I/xz_R1.fq.xz; xz -c $A2 > $I/xz_R2.fq.xz; }
zstd -qc $A1 > $I/zst_R1.fq.zst; zstd -qc $A2 > $I/zst_R2.fq.zst
# one read of 3 MB (longer than the 1 MB inflate block) between normal reads, plain and gz
python3 - <<'EOF' > $I/long.fq
import random
r = random.Random(5)
src = open("/home/falk/audit6/io/reads/sa_R1.fq").read().split("\n")
recs = ["\n".join(src[i:i+4]) for i in range(0, 400, 4)]
L = 3 * 1024 * 1024
s = "".join(r.choice("ACGT") for _ in range(L))
out = recs[:50] + ["@long\n" + s + "\n+\n" + "I" * L] + recs[50:]
print("\n".join(out))
EOF
gzip -c $I/long.fq > $I/long.fq.gz
python3 $S/mkreads.py bgzf $I/long.fq $I/longb.fq.gz

runse() { # runse NAME R1 [extra]
  local name=$1 r1=$2; shift 2
  prun $name.log --db $DB -1 $r1 --prefix $name -o out -t 2 --no_qcmsa --no_strains --sam_format sam --model_se $DB/model_pe.xml "$@" > $name.rc
  local sam=out/$name.sam rec=0 q=0
  if [ -e $sam ]; then records $sam > $name.rec; rec=$(wc -l < $name.rec); q=$(cut -f1 $name.rec | sort -u | wc -l); fi
  printf '%-14s %s sam=%s records=%s qnames=%s | %s\n' $name "$(cat $name.rc)" "$([ -e $sam ] && echo yes || echo no)" $rec $q \
     "$(grep -a -i -E 'error|abort|truncat|corrupt|fail|differ|malformed|unrecognized|warn' $name.log | grep -a -v -E '^(snp|msa)' | head -3 | tr '\n' '|' | cut -c1-300)"
}
runpe() {
  local name=$1 r1=$2 r2=$3; shift 3
  prun $name.log --db $DB -1 $r1 -2 $r2 --prefix $name -o out -t 2 --no_qcmsa --no_strains --sam_format sam "$@" > $name.rc
  local sam=out/$name.sam rec=0 same=-
  if [ -e $sam ]; then records $sam > $name.rec; rec=$(wc -l < $name.rec); same=$(cmp -s $name.rec $BASE && echo SAME || echo diff); fi
  printf '%-14s %s sam=%s records=%s vs_base=%s | %s\n' $name "$(cat $name.rc)" "$([ -e $sam ] && echo yes || echo no)" $rec $same \
     "$(grep -a -i -E 'error|abort|truncat|corrupt|fail|differ|malformed|unrecognized|warn' $name.log | grep -a -v -E '^(snp|msa)' | head -3 | tr '\n' '|' | cut -c1-300)"
}
echo "--- single-end reference (plain R1)"
runse se_base $A1
for v in gz multi bgzf trunc corrupt crc btrunc bnoeof bcrc mcut trail emptygz bgzfplusgz; do runse se_$v $I2/${v}_R1.fq.gz; done
for v in notfq empty blank lower crlf; do runse se_$v $I2/${v}_R1.fq; done
runse se_fasta $I2/fasta_R1.fa
runse se_mbad $I/mbad_R1.fq.gz
[ -e $I/bz_R1.fq.bz2 ] && runse se_bz $I/bz_R1.fq.bz2
[ -e $I/xz_R1.fq.xz ] && runse se_xz $I/xz_R1.fq.xz
runse se_zst $I/zst_R1.fq.zst
runse se_long $I/long.fq
runse se_longgz $I/long.fq.gz
runse se_longb $I/longb.fq.gz
echo "records per variant vs se_base:"
for v in gz multi bgzf crlf trail mbad mcut; do [ -e se_$v.rec ] && echo "  se_$v $(cmp -s se_$v.rec se_base.rec && echo SAME || echo diff)"; done
for v in long longgz longb; do [ -e se_$v.rec ] && echo "  se_$v qnames: $(cut -f1 se_$v.rec | sort -u | wc -l)"; done
echo "--- paired: other compressors, damaged 2nd member, lowercase via the old whole-read path"
runpe mbad $I/mbad_R1.fq.gz $I/mbad_R2.fq.gz
[ -e $I/bz_R1.fq.bz2 ] && runpe bz $I/bz_R1.fq.bz2 $I/bz_R2.fq.bz2
runpe zst $I/zst_R1.fq.zst $I/zst_R2.fq.zst
runpe notfq $I2/notfq_R1.fq $I2/notfq_R2.fq
runpe lower_whole $I2/lower_R1.fq $I2/lower_R2.fq --whole_read_alignment
runpe upper_whole $A1 $A2 --whole_read_alignment
