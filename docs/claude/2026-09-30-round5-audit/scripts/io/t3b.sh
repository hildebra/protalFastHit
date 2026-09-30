#!/bin/bash
# T3 (continued after the WSL restart): the rest of t3_more.sh, one protal at a time.
set -u
uptime
source $(dirname "$0")/lib.sh
T=$W/t3; cd $T
I=$T/in; I2=$W/t2/in
A1=$R/sa_R1.fq; A2=$R/sa_R2.fq
BASE=$W/t1/sam_t1.rec
n=$(wc -l < $A1)
eval "$(sed -n '/^runse() {/,/^}/p; /^runpe() {/,/^}/p' $S/t3_more.sh)"
runse se_xz $I/xz_R1.fq.xz
runse se_zst $I/zst_R1.fq.zst
/usr/bin/time -f "maxrss_kb=%M wall=%e" -o se_long.time true
runse se_long $I/long.fq
runse se_longgz $I/long.fq.gz
runse se_longb $I/longb.fq.gz
echo "records per variant vs se_base:"
for v in gz multi bgzf crlf trail mbad mcut fasta; do [ -e se_$v.rec ] && echo "  se_$v $(cmp -s se_$v.rec se_base.rec && echo SAME || echo diff) $(wc -l < se_$v.rec)"; done
for v in long longgz longb; do [ -e se_$v.rec ] && echo "  se_$v qnames: $(cut -f1 se_$v.rec | sort -u | wc -l) long: $(grep -c '^long' se_$v.rec)"; done
echo "--- paired"
runpe mbad $I/mbad_R1.fq.gz $I/mbad_R2.fq.gz
runpe mbad1 $I/mbad_R1.fq.gz $I2/multi_R2.fq.gz
runpe bz $I/bz_R1.fq.bz2 $I/bz_R2.fq.bz2
runpe zst $I/zst_R1.fq.zst $I/zst_R2.fq.zst
runpe notfq $I2/notfq_R1.fq $I2/notfq_R2.fq
runpe lower_whole $I2/lower_R1.fq $I2/lower_R2.fq --whole_read_alignment
echo "gzip CLI on the damaged multi-member file:"; gzip -t $I/mbad_R1.fq.gz; echo "gzip -t exit $?"; echo "$(gzip -dc $I/mbad_R1.fq.gz 2>/dev/null | wc -l) of $n lines"
