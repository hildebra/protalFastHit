#!/usr/bin/env bash
# Builds protal (with protal_tests) in ~/det-order/<name>/tree/build.
#   build.sh <name> <source> [reorder]
# source: a commit (built from `git archive`) or "tree" (the checkout's working tree, without data/docs/local);
# reorder: also apply scramble_maps.pl, the experiment that makes the maps keyed by taxid iterate in another order.
set -uo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/det-order; NAME=$1; SOURCE=$2; REORDER=${3:-}
E=$W/$NAME; mkdir -p $E
if [ "$SOURCE" = tree ]; then
  rsync -a --checksum --no-times --delete --exclude '/.git' --exclude '/build*/' --exclude '/data' --exclude '/docs' --exclude '/local' $REPO/ $E/tree/
else
  rm -rf $E/src.new && mkdir -p $E/src.new && git -C $REPO archive $SOURCE | tar -x -C $E/src.new
  rsync -a --checksum --no-times --delete --exclude '/build*/' $E/src.new/ $E/tree/ && rm -rf $E/src.new
fi
if [ -n "$REORDER" ]; then perl $HERE/scramble_maps.pl $E/tree || exit 1; fi
cd $E/tree && { [ -d build ] || nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/$NAME.configure.log 2>&1; } && \
  nice cmake --build build --target protal protal_tests -j 6 > $W/$NAME.build.log 2>&1 && echo "OK build $NAME" || { echo "FAIL $NAME"; grep -E 'error' -A3 $W/$NAME.build.log | head -60; exit 1; }
