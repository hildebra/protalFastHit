#!/bin/bash
# Build protal at COMMIT into ~/audit5/builds/<LABEL>; binary and qcmsa to ~/audit5/bin/{protal,qcmsa}_<LABEL>.
# Usage: build_commit.sh COMMIT LABEL
set -eu
C=$1; L=$2
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
D=$HOME/audit5/builds/$L
rm -rf $D; mkdir -p $D
git -C $REPO archive $C | tar -x -C $D
cd $D
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release > cmake.log 2>&1 || { tail -20 cmake.log; exit 1; }
cmake --build build --target protal -j 4 > build.log 2>&1 || { grep -m5 error build.log; exit 1; }
cp build/protal $HOME/audit5/bin/protal_$L
cp scripts/qcmsa.py $HOME/audit5/bin/qcmsa_$L.py
echo "built $C as $L"
