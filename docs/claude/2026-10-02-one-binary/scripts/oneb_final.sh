#!/usr/bin/env bash
# The final change (oneb_final.patch on d11381f): check.sh (clean tree, Release build with tests, ctest, mini
# database, e2e), callgrind of that build (cl3) with the comparison of its outputs, both dispatch paths under gdb,
# and `just install` into a scratch prefix that holds an earlier install's protal_baseline and protal_avx2.
# CG=0 skips the callgrind runs.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
bash $S/check.sh oneb-final $S/oneb_final.patch d11381f || { echo "FAIL check"; exit 1; }
B=$HOME/mt-work/oneb-final/src/build/protal
echo "clone symbols: $(nm $B | grep -c '\.arch_x86_64_v3$') v3, $(nm $B | grep -c '\.default$') default, $(nm $B | grep -c '\.resolver$') resolvers"
CL=$B bash $S/oneb_nocpu.sh 2>/dev/null | grep -vE "^awk|^ *[0-9a-f]+ [0-9a-f]+ [a-zA-Z] __cpu"
CL=$B bash $S/oneb_v3cpu.sh | grep -E "already hit|same|DIFFER" | tr '\n' ' '; echo
W=$HOME/mt-work/onebcg; I=$HOME/mt-work/isacg; DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
declare -A ARGS=(
  [pe100k]="-1 $I/r1.fq -2 $I/r2.fq --no_profile"
  [ont3M]="-1 $P/ont_b3000000/sim/reads/ont_b3000000_s_1.fq.gz --read_type ont --no_profile"
  [prof500k]="--profile_only $HOME/bench071/runs/v071.full.pe.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst" )
one() {
  local w=$1 b=cl3; rm -rf $W/out.$w.$b
  nice valgrind --tool=callgrind --callgrind-out-file=$W/cg.$w.$b $B --db $DB ${ARGS[$w]} --prefix s -o $W/out.$w.$b -t 1 --no_qcmsa > $W/log.$w.$b 2>&1
  callgrind_annotate --threshold=100 $W/cg.$w.$b 2>/dev/null | c++filt > $W/$w.$b.txt
  echo "done $w $b: $(grep 'PROGRAM TOTALS' $W/$w.$b.txt)"
  r=$(diff -r -q -x '*_runtime.tsv' -x '*.sam.zst' $W/out.$w.base $W/out.$w.$b > /dev/null && echo same || echo DIFFER)
  s=""; [ -f $W/out.$w.base/s.sam.zst ] && { cmp -s <(zstdcat $W/out.$w.base/s.sam.zst) <(zstdcat $W/out.$w.$b/s.sam.zst) && s="SAM text same" || s="SAM text DIFFER"; }
  echo "$w $b vs base: outputs $r $s"
}
[ "${CG:-1}" = 1 ] && { for w in pe100k ont3M prof500k; do one $w & done; wait; }
# just install over an earlier install's layout
PF=$HOME/mt-work/oneb-install; rm -rf $PF; mkdir -p $PF/bin
for f in protal protal_baseline protal_avx2; do printf '#!/bin/sh\necho old %s\n' $f > $PF/bin/$f; chmod +x $PF/bin/$f; done
(cd $HOME/mt-work/oneb-final/src && nice just install $PF > $HOME/mt-work/oneb-final/install.log 2>&1) || { echo "FAIL just install"; tail -5 $HOME/mt-work/oneb-final/install.log; }
tail -1 $HOME/mt-work/oneb-final/install.log
echo "installed: $(ls $PF/bin | tr '\n' ' ')"; file -b $PF/bin/protal | cut -c1-60
echo FINAL DONE
