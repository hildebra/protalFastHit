#!/bin/bash
# Compiles and runs chunk_cap.cpp against the (unmodified) strain-build headers, with the include
# directories and definitions of the protal target.
S=$(dirname "$0")
W=~/audit6/longreads
SRC=~/strain-build/src
mkdir -p $W/chunkcap && cd $W/chunkcap
cp $S/chunk_cap.cpp .
F=$W/build-asan/CMakeFiles/protal.dir/flags.make
INC=$(grep CXX_INCLUDES $F | sed "s/CXX_INCLUDES = //" | sed "s#$W/src#$SRC#g")
DEF=$(grep CXX_DEFINES $F | sed "s/CXX_DEFINES = //")
g++ -std=c++20 -O1 -fopenmp $DEF $INC -include iostream chunk_cap.cpp -o chunk_cap > build.log 2>&1 || { grep -m5 error build.log; exit 1; }
./chunk_cap > $S/out_chunkcap.txt 2>&1
echo "rc=$?"
