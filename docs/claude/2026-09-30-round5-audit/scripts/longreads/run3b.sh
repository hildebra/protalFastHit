#!/bin/bash
# The four reads again with --x_drop 0 (WFA2's X-drop off; wf-adaptive stays on), release binary.
S=$(dirname "$0")
W=~/audit6/longreads
P=~/strain-build/bin/protal
PA=$W/build-asan/protal
D=$W/mini_db
OUT=$S/out_e3b.txt
cd $W/e2c
{
for xd in 1000 0; do
  $P --db $D -1 four.fq --read_type ont --model_ont $D/model_pe.xml -o out_xd$xd -t 1 --sam_format sam --no_qcmsa --no_strains --no_profile -m 10 --prefix four --x_drop $xd > xd$xd.log 2>&1
  echo "--x_drop $xd rc=$? ; records on genes 35/38/103/82 at read ends:"
  grep -v '^@' out_xd$xd/four.sam | awk '($3 ~ /_(35|38|103|82)$/) {h=$6; gsub(/[0-9]+[MXID]/,"",h); print "   ", $1,$2,$3,$4,$5,length($10),h}'
done
echo "== debug build, --x_drop 0"
PROTAL_LR_DEBUG=1 ASAN_OPTIONS=detect_leaks=0 $PA --db $D -1 four.fq --read_type ont --model_ont $D/model_pe.xml -o outdbg0 -t 1 --sam_format sam --no_qcmsa --no_strains --no_profile -m 10 --prefix four --x_drop 0 > four_dbg0.log 2>&1
echo "rc=$?"
grep -E "^LRDBG (read|align-fail)" four_dbg0.log | cut -c1-300
grep -E "^LRDBG (read|segment)" four_dbg0.log | grep -E "read|_35:|_38:|_103:|_82:" | cut -c1-300
grep -E "^(Invalid alignment)" four_dbg0.log | cut -c1-300
} > $OUT 2>&1
echo done
