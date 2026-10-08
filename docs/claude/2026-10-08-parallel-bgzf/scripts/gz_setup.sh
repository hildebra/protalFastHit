#!/usr/bin/env bash
# Experiment set-up: HEAD's source in ~/gzpar/src, readbench built against it, the 5M-pair sample as BGZF (as is),
# one gzip member (gzip -6) and uncompressed. 4 cores, niced.
set -uo pipefail
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
G=$HOME/gzpar; T="taskset -c 0-3 nice -n 5"
B=$HOME/bench071/samples_deep/points/rl150_p5000000/sim/reads
mkdir -p $G/data $G/rb
rm -rf $G/src; mkdir -p $G/src
git -C $REPO archive HEAD | tar -x -C $G/src
git -C $REPO rev-parse --short HEAD > $G/src.commit
cp $REPO/docs/claude/2026-10-01-multithreading-audit/scripts/readbench.cpp $G/rb/
cd $G/rb
S=$G/src
$T g++ -std=c++20 -O3 -march=x86-64-v2 -fopenmp -I$S/src -I$S/src/IO -I$S/src/SequenceUtils -I$S/src/Utilities -I$S/lib \
   readbench.cpp $S/src/IO/FastxReader.cpp -lisal -lzstd -lpthread -o readbench 2> build.err || { head -40 build.err; exit 1; }
echo "readbench built ($(cat $G/src.commit))"
for r in R1 R2; do
  [ -s $G/data/$r.fq ] || $T bash -c "zcat $B/rl150_p5000000_s_1_$r.fq.gz > $G/data/$r.fq"
  [ -s $G/data/$r.member.gz ] || $T bash -c "gzip -6 -c $G/data/$r.fq > $G/data/$r.member.gz"
  ln -sf $B/rl150_p5000000_s_1_$r.fq.gz $G/data/$r.bgzf.gz
done
ls -laL $G/data
for f in $G/data/R1.bgzf.gz $G/data/R1.member.gz; do echo "$f: $(head -c 18 $f | od -An -tx1 | tr -s ' ')"; done
echo SETUP_DONE
