#!/bin/bash
# protal and simulate_metagenomes of one commit, in a clone of the Windows checkout on WSL's file system (so that the
# binaries' --version names the commit and build_gtdb_database.py's check finds its scripts at the same commit).
# usage: setup_tools.sh <commit> <folder>     (WSL; 4 cores)
set -eu
COMMIT=$1; DIR=$2
if [ ! -d "$DIR/.git" ]; then
  git clone --quiet /mnt/c/Users/hildebra/Documents/locDev/protal "$DIR"
fi
cd "$DIR"
git fetch --quiet origin
git checkout --quiet "$COMMIT"
echo "source at $(git log -1 --format='%h %s' | cut -c1-100)"
# The system's compilers and CMake (the protal-db-build environment's CMake does not look in /usr for ISA-L).
[ -x build/protal ] || rm -rf build
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > build_configure.log
taskset -c 0-3 nice -n 5 cmake --build build -j 4 --target protal simulate_metagenomes > build.log
./build/protal --version | tail -1
./build/simulate_metagenomes --version | tail -1
