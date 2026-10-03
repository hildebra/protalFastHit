#!/usr/bin/env bash
# Sync the denoise worktree into WSL and build protal, simulate_metagenomes and the unit tests (Release).
set -euo pipefail
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal-denoise
B=$HOME/protal-denoise
mkdir -p $B/src
rsync -a --checksum --no-times --delete --exclude /.git --exclude '/build*/' --exclude /data --exclude /local \
  --exclude '__pycache__' $SRC/ $B/src/
cd $B
if [ ! -f build/build.ninja ]; then
  cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > build_configure.log 2>&1 \
    || { tail -30 build_configure.log; exit 1; }
fi
nice -n 5 cmake --build build --target ${TARGETS:-protal simulate_metagenomes protal_tests} -j ${JOBS:-6} > build.log 2>&1 \
  || { grep -E "error|Error" build.log | head -60; tail -20 build.log; exit 1; }
tail -3 build.log
ls -la build/protal build/simulate_metagenomes 2>/dev/null || find build -maxdepth 2 -name protal -type f
