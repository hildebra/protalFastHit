#!/usr/bin/env bash
# replay.sh [ROUNDS] [VARIANTS...]: the dumped calls described, then wfareplay (built against ~/wfasimd/src) on each
# sample with the variants alternated; the per-call results compared between variants.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfasimd
W=$HOME/wfasimd; H=${H:-$W/src}; C=$W/calls; R=$W/replay; BIN=${BIN:-$R/wfareplay}
ROUNDS=${1:-7}; shift || true
VARIANTS=${*:-classic singletrack classic/2 singletrack/2}
mkdir -p $R
for s in pe pb ont; do
  awk -F'\t' -v s=$s '{n[$1]++; q+=length($8); r+=length($9); if ($7 < 2147483647) m+=$7} END {for (k in n) printf "%s: %s %d calls; ", s, k, n[k]; printf "mean read part %.0f, reference %.0f, budget %.0f\n", q/NR, r/NR, m/NR}' $C/$s.tsv
done
if [ ! -x $BIN ] || [ $S/wfareplay.cpp -nt $BIN ]; then
  g++ -O3 -march=x86-64 -mtune=generic -std=c++20 -I$H/src -I$H/src/IO -I$H/src/Hash -I$H/src/Core -I$H/src/SequenceUtils -I$H/src/Alignment \
    -I$H/src/Utilities -I$H/src/Taxonomy -I$H/lib -I$H/lib/robin -I$H/lib/tsl -I$H/lib/wfa2-lib -I$H/lib/gzstream -I$H/build/zlib-ng \
    -I$H/lib/zlib-ng -I$H/build/generated $S/wfareplay.cpp $H/build/src/libwfa_lib.a -o $BIN 2>&1 | grep -E ' error' | head
fi
for s in pe pb ont; do
  echo "-- $s"
  taskset -c 0-3 $BIN $C/$s.tsv $ROUNDS $R/$s $VARIANTS 2>/dev/null | sed 's/^/   /'
done
for s in pe pb ont; do
  b=$R/$s.classic; [ -f $b ] || continue
  for v in classic2 singletrack singletrack2; do
    f=$R/$s.$v; [ -f $f ] || continue
    same=$(paste -d' ' $b $f | awk '$1==$4 && $2==$5' | wc -l)
    ops=$(paste -d' ' $b $f | awk '$1==$4 && $2==$5 && $3!=$6' | wc -l)
    echo "   $s classic vs $v: outcome and score the same in $same of $(wc -l < $b), other CIGAR in $ops"
  done
  [ -f $R/$s.singletrack2 ] && echo "   $s singletrack vs singletrack/2: $(cmp -s $R/$s.singletrack $R/$s.singletrack2 && echo identical || echo DIFFER)"
done
