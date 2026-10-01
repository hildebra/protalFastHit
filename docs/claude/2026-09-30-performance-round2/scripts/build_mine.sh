#!/bin/bash
# A build of HEAD with only the listed working-tree files on top (other work in the same checkout stays out):
#   build_mine.sh FILE...   (paths relative to the checkout) -> $PERF_DIR/pair/mine/build
# The first run unpacks HEAD; later runs copy the files again (cp sets their mtime to now, so ninja rebuilds them;
# rsync -a would keep the old mtimes, and a build could look up to date). FRESH=1 starts from HEAD again.
set -e
source "$(dirname "$0")/env.sh"
M=$PERF_DIR/pair/mine
if [ -n "${FRESH:-}" ] || [ ! -d $M/src ]; then
  rm -rf $M; mkdir -p $M/src
  (cd $PROTAL_SRC && git archive HEAD) | tar -x -C $M/src
fi
for f in "$@"; do mkdir -p $M/src/$(dirname $f); cp $PROTAL_SRC/$f $M/src/$f; done
cd $M
[ -f build/build.ninja ] || cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > cmake.log 2>&1
cmake --build build --target protal protal_avx2 simulate_metagenomes protal_tests -j ${JOBS:-4} > build.log 2>&1 && echo "mine built" || { grep -m5 "error:" build.log; exit 1; }
