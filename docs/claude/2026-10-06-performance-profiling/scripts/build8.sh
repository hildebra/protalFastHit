#!/usr/bin/env bash
# build8.sh ref|work|isal: ~/perf8/ref = git archive HEAD; ~/perf8/work = the same with this session's changed files from
# the working tree copied over it (only the files listed in ~/perf8/files.txt, so the other sessions' uncommitted
# edits stay out); ~/perf8/isal = work plus the ISA-L files (~/perf8/files_isal.txt), built against the ISA-L in
# ~/isal/prefix (no libisal-dev here). Release with the unit tests.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/perf8; N=$1; E=$W/$N
mkdir -p $W
if [ "$N" = ref ] || [ ! -d $E/build ]; then
  rm -rf $E; mkdir -p $E
  git -C $REPO archive HEAD | tar -x -C $E
  git -C $REPO rev-parse --short HEAD > $W/$N.commit
fi
lists=""; extra=""
[ "$N" = work ] && lists="$W/files.txt"
[ "$N" = isal ] && { lists="$W/files.txt $W/files_isal.txt"; extra="-DCMAKE_PREFIX_PATH=$HOME/isal/prefix"; }
for l in $lists; do
  while read -r f; do [ -n "$f" ] && cp "$REPO/$f" "$E/$f" && touch "$E/$f"; done < $l
done
cd $E
[ -d build ] || nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON $extra > $W/$N.configure.log 2>&1
nice cmake --build build -j 5 > $W/$N.build.log 2>&1 && echo "OK build $N ($(cat $W/$N.commit))" || { echo "FAIL build $N"; grep -E 'error|Error' -A4 $W/$N.build.log | head -60; }
