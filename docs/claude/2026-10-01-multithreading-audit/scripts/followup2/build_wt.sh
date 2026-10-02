#!/usr/bin/env bash
# build_wt.sh NAME PATCH [BASE] [TARGETS]: a64b60f-based tree + PATCH in ~/mt-work/NAME, build TARGETS (default protal).
set -uo pipefail
NAME=$1; PATCH=$2; BASE=${3:-a32b60f}; TARGETS=${4:-protal}
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/mt-work/$NAME
rm -rf $W/src.new; mkdir -p $W/src.new
git -C $REPO archive $BASE | tar -x -C $W/src.new
(cd $W/src.new && patch -p1 -s < $PATCH) || { echo "FAIL patch"; exit 1; }
mkdir -p $W/src
rsync -a --checksum --no-times --delete --exclude /build $W/src.new/ $W/src/
rm -rf $W/src.new
cd $W/src
[ -f build/build.ninja ] || nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/configure.log 2>&1
nice cmake --build build --target $TARGETS -j 4 > $W/build.log 2>&1 || { echo "FAIL build"; grep -E 'error' -A3 $W/build.log | head -60; exit 1; }
echo "OK build $NAME"
