#!/bin/bash
# protal at the working tree (~/pgzfix/build/protal) on the mini DB (~/pgzfix/mini_db, scripts/mini_db/build_mini_db.sh)
# with the same reads as one gzip member each and as members of 64 KiB, -t 16 (2 inflating threads per file),
# --verbose: how each file is read, and whether the outputs are the same.
set -uo pipefail
SP=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-gzip-members-regression/scripts
W=~/pgzfix
C=$W/check
rm -rf $C
mkdir -p $C
cd $C
zstd -dc $W/mini_db/protal_db/full_reference.fna.zst > ref.fna
perl $SP/mini_reads.pl ref.fna 400000 3 reads
for m in 1 2; do
  gzip -6 -c reads_R$m.fq > one_R$m.fq.gz
  split -b 64K --filter='gzip -6' reads_R$m.fq > members_R$m.fq.gz
done
ls -la *.gz
for kind in one members; do
  taskset -c 0-3 nice -n 5 $W/build/protal --db $W/mini_db/protal_db -t 16 --verbose -1 ${kind}_R1.fq.gz -2 ${kind}_R2.fq.gz \
    --prefix $kind -o out_$kind > $kind.log 2>&1
  echo "$kind exit $?"
  grep -E "^Input |Aligning reads took|records written" $kind.log
done
for f in $(cd out_one && find . -name '*.profile*' -type f | sort); do
  g=${f//one/members}
  if cmp -s out_one/$f out_members/$g; then echo "same: $f"; else echo "DIFFERENT: $f $g"; fi
done
