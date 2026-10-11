#!/usr/bin/env bash
# issue_repro.sh: wfa_issue_repro.c against three WFA2-lib trees: upstream v2.3.6, the fork tmokveld/WFA2-lib 85d063f
# as it is, and protal's patched copy of the fork (lib/wfa2-lib of the merge being checked). Each library from its C
# sources with gcc -O2, as protal builds it (no -march=native).
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad
F=$S/dl2/fork; W=$HOME/issuerepro
rm -rf $W; mkdir -p $W/upstream $W/fork
git -C $F archive v2.3.6 | tar -x -C $W/upstream
git -C $F archive 85d063f | tar -x -C $W/fork
cp -r $HOME/mcheck/src/lib/wfa2-lib $W/patched
SRC=$(sed -n '/set(WFA2_SOURCES/,/)/p' $HOME/mcheck/src/lib/wfa2-lib.cmake | grep -oE '[a-z_0-9]+/[a-z_0-9]+\.c')
for t in upstream fork patched; do
  cd $W/$t && mkdir -p obj
  for f in $SRC; do gcc -O2 -I. -c $f -o obj/$(basename $f .c).o 2> /dev/null || echo "compile failed: $t $f"; done
  ar rcs libwfa.a obj/*.o
  def=""; [ $t != upstream ] && def="-DHAVE_SINGLETRACK"
  gcc -O2 $def -I. $S/wfa_issue_repro.c libwfa.a -lm -o repro 2>&1 | grep -E "error" | head -3
done
for c in "1 high" "1 singletrack" "2 high" "2 singletrack"; do
  for t in upstream fork patched; do
    echo "-- $t: $c"
    (cd $W/$t && ./repro $c 2>&1 | grep -v '^[<>]'; echo "   exit ${PIPESTATUS[0]}")
  done
done
