#!/bin/bash
# The SAM records (sorted, without the header) of the base and new runs, round 1: their count and hash.
cd ~/keysfix/perf
for t in pe hifi; do for keys in direct compact; do for b in base new; do
  f=seed_${b}_${t}_${keys}_1/${t}.sam.zst
  printf "%s\t%s\n" "$f" "$(zstd -dcq $f | grep -v '^@' | sort | tee >(wc -l > /tmp/keysfix_n) | md5sum | cut -c1-12) $(sleep 0.2; cat /tmp/keysfix_n)"
done; done; done
for c in sparse sparsec; do for b in base new; do
  f=load_${b}_${c}_1/s.sam.zst
  printf "%s\t%s\n" "$f" "$(zstd -dcq $f | grep -v '^@' | sort | md5sum | cut -c1-12)"
done; done
rm -f /tmp/keysfix_n
