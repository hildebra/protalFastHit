#!/usr/bin/env bash
# Build and run the ReadFrame allocation probe. Plain build under a ~1.2 GB address-space cap:
# a 4 GB resize throws std::bad_alloc, proving the pre-validation allocation.
set -u
S=$(dirname "$0")
W=~/audit6/database
mkdir -p $W/hbuild
cp $S/harness.cpp $W/hbuild/
cp /mnt/c/Users/hildebra/Documents/locDev/protal/.claude/worktrees/strain-fixes/src/Utilities/Zstd.h $W/hbuild/
cd $W/hbuild
g++ -O2 -std=c++20 -I. harness.cpp -o harness -lzstd 2>err.log
if [ ! -x harness ]; then echo "BUILD FAILED"; cat err.log; exit 1; fi
echo "=== run under 1.2 GB address-space cap (ulimit -v)"
/usr/bin/env bash -c "ulimit -v 1200000; ./harness" 2>&1
