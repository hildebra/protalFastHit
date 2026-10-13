#!/bin/bash
# Test files for the before/after benchmark: R1.fq (make_fastq.pl 4000000 7) as one member, members of 8 MiB, 1 MiB and
# 64 KiB of FASTQ, and mixed.gz: its first 400 MB as one member, the rest in members of 64 KiB.
set -uo pipefail
SP=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-gzip-members-regression/scripts
D=~/pgzbench/data
mkdir -p $D
cd $D
[ -s R1.fq ] || taskset -c 0-3 nice -n 5 perl $SP/make_fastq.pl 4000000 7 > R1.fq
[ -s one.gz ] || taskset -c 0-3 nice -n 5 gzip -6 -c R1.fq > one.gz
[ -s m64k.gz ] || taskset -c 0-3 nice -n 5 split -b 64K --filter='gzip -6' R1.fq > m64k.gz
[ -s m1m.gz ] || taskset -c 0-3 nice -n 5 split -b 1M --filter='gzip -6' R1.fq > m1m.gz
[ -s m8m.gz ] || taskset -c 0-3 nice -n 5 split -b 8M --filter='gzip -6' R1.fq > m8m.gz
if [ ! -s mixed.gz ]; then
  head -c 400000000 R1.fq | taskset -c 0-3 nice -n 5 gzip -6 > mixed.gz
  tail -c +400000001 R1.fq | taskset -c 0-3 nice -n 5 split -b 64K --filter='gzip -6' >> mixed.gz
fi
ls -la
