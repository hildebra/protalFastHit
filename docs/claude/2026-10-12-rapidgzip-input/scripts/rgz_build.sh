#!/bin/bash
# Build protal with lib/isa-l + lib/rapidgzip in ~/rgzbuild (copy of the working tree), 4 cores; run the unit tests.
set -uo pipefail
P=/mnt/c/Users/hildebra/Documents/locDev/protal
W=~/rgzbuild
mkdir -p $W/src
rsync -a --checksum --no-times --delete \
  --exclude='/data' --exclude='/build*/' --exclude='/local' --exclude='/.git' --exclude='/.claude' --exclude='/misc' \
  --exclude='/tests/misc' --exclude='/docs/claude' \
  $P/ $W/src/
cd $W
if [ ! -f build/build.ninja ]; then
  taskset -c 0-3 nice -n 5 cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > configure.log 2>&1
  echo "configure exit $?"
  grep -i "isa-l\|nasm\|rapidgzip\|error" configure.log | head -20
fi
taskset -c 0-3 nice -n 5 cmake --build build --target protal simulate_metagenomes protal_tests -- -j4 > build.log 2>&1
echo "build exit $?"
grep -E "error|Error|warning: .*rapidgzip" build.log | head -40
tail -3 build.log
