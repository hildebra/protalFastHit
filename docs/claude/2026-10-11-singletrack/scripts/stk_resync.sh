#!/usr/bin/env bash
# stk_resync.sh: ~/stk/work as 7d5df4b + st2.patch (files the patch touches restored from the archive first), rebuilt.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/stk/work
files=$(grep -E '^diff --git' $S/st2.patch | awk '{print substr($3,3)}' | sort -u)
fresh=$HOME/stk/fresh; rm -rf $fresh; mkdir -p $fresh; git -C $REPO archive 7d5df4b | tar -x -C $fresh
for f in $files; do
  if [ -e $fresh/$f ]; then cp $fresh/$f $W/$f; else rm -f $W/$f; fi
done
# lib/wfa2-lib: start from v2.3.6 as in the archive (the patch replaces it)
rm -rf $W/lib/wfa2-lib && cp -r $fresh/lib/wfa2-lib $W/lib/wfa2-lib
(cd $W && git apply --whitespace=nowarn $S/st2.patch) || { echo "FAIL patch"; exit 1; }
cp $S/st2.patch $S/st.patch
rm -rf $fresh
(taskset -c 0-3 nice -n 5 cmake --build $W/build -j 4 > $HOME/stk/work.build3.log 2>&1) && echo "OK build work" \
  || { echo "FAIL build work"; grep -E 'error' -A4 $HOME/stk/work.build3.log | head -40; }
echo "BUILDDONE $(date +%T)"
