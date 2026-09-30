#!/usr/bin/env bash
# Build (ASan) and run the DecodeChunk fuzz. Needs IndexCodec.h + Zstd.h include path.
set -u
SRC=/mnt/c/Users/hildebra/Documents/locDev/protal/.claude/worktrees/strain-fixes/src
S=$(dirname "$0")
W=~/audit6/database
mkdir -p $W/fuzz
# Copy the two headers into a private tree so nothing on /mnt/c is built in place.
mkdir -p $W/fuzz/inc/Hash $W/fuzz/inc/Utilities
cp $SRC/Hash/IndexCodec.h $W/fuzz/inc/Hash/
cp $SRC/Utilities/Zstd.h $W/fuzz/inc/Utilities/
# IndexCodec.h includes "Zstd.h" (same dir); provide it next to it too.
cp $SRC/Utilities/Zstd.h $W/fuzz/inc/Hash/
cp $S/fuzz_decode.cpp $W/fuzz/
cd $W/fuzz
g++ -O1 -g -std=c++20 -fsanitize=address,undefined -fno-omit-frame-pointer -I inc fuzz_decode.cpp -o fuzz -lzstd 2>build.log
if [ ! -x fuzz ]; then echo BUILD FAILED; tail -30 build.log; exit 1; fi
echo "=== run"
ASAN_OPTIONS=abort_on_error=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 ./fuzz 2>&1 | tail -25
echo "fuzz exit=${PIPESTATUS[0]}"
