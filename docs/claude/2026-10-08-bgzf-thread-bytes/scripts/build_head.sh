#!/bin/bash
# A clean WSL tree of 401c4f5 (head.tar), configured as ~/pipefix (Release, Ninja, ccache, tests), and its unit test
# binary built on cores 0-3; the binary is kept as ~/bgzfdet/bin/tests_head.
set -u
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/99a2c9fe-d722-41e1-bc32-3fc6baabdeb2/scratchpad
W=~/bgzfdet
mkdir -p $W/src $W/bin
cd $W/src
tar -xf $SP/head.tar
echo "extracted $(date +%T)"
cd $W
taskset -c 0-3 nice -n 5 cmake -S src -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
    -DPROTAL_BUILD_TESTS=ON > cmake.log 2>&1
echo "cmake rc $? $(date +%T)"
taskset -c 0-3 nice -n 5 ninja -C build -j4 protal_tests > ninja_head.log 2>&1
rc=$?
echo "ninja rc $rc $(date +%T)"
tail -3 ninja_head.log
[ $rc -eq 0 ] && cp build/tests/protal_tests bin/tests_head && echo "kept bin/tests_head"
