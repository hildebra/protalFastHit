#!/usr/bin/env bash
# The working tree without the unmapped fast path (its call commented out), built in ~/perf-gtdb/work_nofast.
set -uo pipefail
W=$HOME/perf-gtdb; E=$W/work_nofast
rsync -a --delete --exclude '/build*/' $W/work/ $E/
sed -i 's|^\(\s*\)if (SkipUnmapped(m_tokens)) continue;|\1// no fast path|' $E/src/IO/SamHandler.h
grep -c "no fast path" $E/src/IO/SamHandler.h
cd $E && { [ -d build ] || nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > $W/nofast.configure.log 2>&1; } && nice cmake --build build --target protal -j 6 > $W/nofast.build.log 2>&1 && echo "OK nofast build" || { echo FAIL; tail -30 $W/nofast.build.log; }
cp $E/build/protal $W/samread/protal.nofast
