#!/usr/bin/env bash
# The experiment tree: HEAD of audit-fixes (git archive) in ~/mt-work/perf3/ref, built (protal); a verbose
# one-thread run on the Nanopore 3 Mb and PacBio 3 Mb samples for the anchored / whole-window counts.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
W=$HOME/mt-work/perf3; mkdir -p $W
rm -rf $W/ref; mkdir -p $W/ref
git -C $REPO archive HEAD | tar -x -C $W/ref
git -C $REPO rev-parse --short HEAD > $W/ref/COMMIT
cd $W/ref && nice cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > $W/ref.configure.log 2>&1 && \
  nice cmake --build build --target protal -j 4 > $W/ref.build.log 2>&1 || { echo FAIL build; exit 1; }
DB=$HOME/bench071/V071/protal_db; P=$HOME/bench071/samples/points
for s in ont_b3000000 pb_b3000000; do
  t=ont; [ $s = pb_b3000000 ] && t=pb
  rm -rf $W/v.$s
  $W/ref/build/protal --db $DB -1 $P/$s/sim/reads/${s}_s_1.fq.gz --read_type $t --no_profile --prefix s -o $W/v.$s -t 1 --no_qcmsa --verbose > $W/v.$s.log 2>&1
  echo "== $s"; grep -E "Anchors aligned|Aligning reads took|segments|reads" $W/v.$s.log | head -12
done
