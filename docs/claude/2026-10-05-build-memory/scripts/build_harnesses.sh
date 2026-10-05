#!/usr/bin/env bash
# Compiles the harnesses against one source tree: build_harnesses.sh TREE OUT_DIR
# TREE is a protal checkout on the Linux file system (git archive of a commit, see setup.sh), so that
# <zstd.h> is the system header, not protal's Zstd.h (/mnt/c is case-insensitive).
set -e
tree=$1; out=$2
here=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$out"
inc="-I$tree/src/Utilities -I$tree/src/Hash -I$tree/src/SequenceUtils -I$here"
g++ -std=c++20 -O3 -march=x86-64 -fopenmp $inc "$here/codec_memory.cpp" -o "$out/codec_memory" -lzstd -lpthread
g++ -std=c++20 -O3 -march=x86-64 -fopenmp $inc "$here/suspect_memory.cpp" -o "$out/suspect_memory" -lpthread
echo "built $out/codec_memory $out/suspect_memory from $tree"
