#!/usr/bin/env bash
# Builds this change in WSL from a Linux-side copy: ~/buildmem/new is HEAD (git archive) with the change's
# files laid over it, so that other sessions' uncommitted work in the checkout stays out. Also builds the unit
# tests. setup.sh builds the baseline (~/buildmem/base).
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal
WORK=$HOME/buildmem
FILES="src/Build.h src/RunProtal.h src/SequenceUtils/GenomeLoader.h src/Hash/IndexCodec.h src/Utilities/Zstd.h"
FILES="$FILES ${EXTRA_FILES:-}"
mkdir -p $WORK && cd $WORK || exit 1
if [ ! -d new/src ]; then mkdir -p new && (cd $SRC && git archive --format=tar HEAD) | tar -x -C new || exit 1; fi
for f in $FILES; do cp "$SRC/$f" "new/$f" || exit 1; touch "new/$f"; done
cmake -S new -B new/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > new-cmake.log 2>&1 || { tail -20 new-cmake.log; exit 1; }
nice cmake --build new/build --target protal protal_tests -j 6 > new-build.log 2>&1 || { grep -E "error" -A4 new-build.log | head -60; echo BUILD FAILED; exit 1; }
ls -la new/build/protal
