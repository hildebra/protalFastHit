#!/usr/bin/env bash
# Where zstd comes from, and what the mini DB bundle holds.
set -u
ls ~/strain-build/src/lib | head -40
ls /usr/include/zstd.h 2>&1; ls /usr/lib/x86_64-linux-gnu/libzstd* 2>&1
grep -n -i zstd ~/strain-build/src/CMakeLists.txt | head -20
ls ~/strain-build/src/build 2>&1 | head
ls ~/strain-build/src/build/lib 2>&1 | head
find ~/strain-build/src/build -name 'libzstd*' 2>/dev/null | head
~/strain-build/bin/protal --help 2>&1 | grep -n -E 'db|compress|unpack|model|read_type' | head -40
