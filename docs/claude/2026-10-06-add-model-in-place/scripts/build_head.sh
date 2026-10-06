#!/bin/bash
# Build HEAD (15006c2, without the in-place --add_model) into ~/protal-addmodel/head from
#   git archive --format=tar -o head.tar 15006c2   (run in the checkout; HEAD_TAR=path/to/head.tar)
set -e
W=~/protal-addmodel/head
rm -rf $W && mkdir -p $W/src
tar -xf $HEAD_TAR -C $W/src
cd $W
cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1 || { tail -30 cmake.log; exit 1; }
nice -n 10 cmake --build build -j 6 --target protal > build.log 2>&1 || { grep error build.log | head; exit 1; }
ls -la build/protal
