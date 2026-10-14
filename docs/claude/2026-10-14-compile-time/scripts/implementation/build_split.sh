#!/bin/bash
# Syncs the split-units worktree to ~/split-build/src and builds TARGETS (default: protal) there, on cores 0-3.
# Changed files get the copy's time (--checksum --no-times), so that ninja sees them.
set -o pipefail
W=/mnt/c/Users/hildebra/Documents/locDev/protal-split
B=~/split-build
mkdir -p $B/src
rsync -a --checksum --no-times --delete --exclude=/.git --exclude=/data --exclude='/build*/' --exclude=/misc $W/ $B/src/
if [ ! -f $B/build/build.ninja ]; then
    cmake -S $B/src -B $B/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON -DCMAKE_EXPORT_COMPILE_COMMANDS=ON > $B/configure.log 2>&1 || { tail -30 $B/configure.log; exit 1; }
fi
cd $B/build
start=$(date +%s.%N)
taskset -c 0-3 nice -n 5 ninja -j4 -k 0 ${@:-protal} > $B/build.log 2>&1
status=$?
end=$(date +%s.%N)
echo "ninja exit $status, $(echo "$end - $start" | bc) s"
grep -E 'error|Error' $B/build.log | grep -v '^\[' | head -${ERRLINES:-60}
tail -3 $B/build.log
