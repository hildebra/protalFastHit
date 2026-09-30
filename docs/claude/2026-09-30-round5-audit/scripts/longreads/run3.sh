#!/bin/bash
# Experiment 3: the four reads whose gene at a read end went to another species, through the ASan
# build with PROTAL_LR_DEBUG diagnostics.
S=$(dirname "$0")
W=~/audit6/longreads
PA=$W/build-asan/protal
D=$W/mini_db
OUT=$S/out_e3.txt
cd $W/e2c
ls -la $PA
{
PROTAL_LR_DEBUG=1 ASAN_OPTIONS=detect_leaks=0 $PA --db $D -1 four.fq --read_type ont --model_ont $D/model_pe.xml -o outdbg -t 1 --sam_format sam --no_qcmsa --no_strains -m 10 --prefix four > four_dbg.log 2>&1
echo "rc=$?"
grep -E "^LRDBG (read|align-fail|drop)" four_dbg.log | cut -c1-400
echo "== candidates of the genes in question"
grep -E "^LRDBG (read|cand) " four_dbg.log | grep -E "read|_35 |_38 |_103 |_82 " | cut -c1-300
echo "== segments with them"
grep -E "^LRDBG (read|segment)" four_dbg.log | grep -E "read|_35:|_38:|_103:|_82:" | cut -c1-300
echo "== invalid"
grep -E "^(CIGAR|Invalid alignment)" four_dbg.log | cut -c1-600
grep -E "runtime error|AddressSanitizer" four_dbg.log | head
} > $OUT 2>&1
echo done
