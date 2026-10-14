#!/bin/bash
# Compile-time study: a build tree of HEAD (git archive) without ccache.
# Usage (in WSL): git archive -o head.tar <commit> (Git Bash), then bash setup.sh /mnt/c/.../head.tar
set -e
TAR=${1:?usage: setup.sh <git archive tar of the commit>}
W=~/ctime
rm -rf $W
mkdir -p $W/src
cd $W/src
tar -xf "$TAR"
cmake -S . -B $W/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXPORT_COMPILE_COMMANDS=ON > $W/configure.log 2>&1
echo configured
# Everything main.cpp.o depends on, built first (libraries, generated headers), timed per target.
cd $W/build
/usr/bin/time -v taskset -c 0-3 nice -n 5 ninja -j4 protal_lib cPMML wfa_lib gzstream_lib zlib-ng protal_commit > $W/deps.log 2> $W/deps.time
tail -3 $W/deps.log
grep -E 'Elapsed|Maximum resident' $W/deps.time
grep -o '"command": "[^"]*main.cpp[^"]*"' compile_commands.json | head -2
