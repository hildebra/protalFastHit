#!/bin/bash
# Bases asked for and made, reads and their median length, of each long-read sample of collections.
# usage: long_read_stats.sh COLLECTION_DIR...
for dir in "$@"; do
  for samples in "$dir"/points/{pb,ont}_b*/sim/samples.tsv; do
    [ -f "$samples" ] || continue
    point=$(basename "$(dirname "$(dirname "$samples")")")
    tail -n +2 "$samples" | while IFS=$'\t' read -r sample reads truth community; do
      zcat "$reads" | awk 'NR % 4 == 2 { print length($0) }' | sort -n | awk -v s="$sample" -v asked="${point#*_b}" \
        '{ a[NR] = $1; b += $1 } END { printf "%-26s asked %11d, got %11d bases (%.2fx) in %6d reads, median %6d bp\n", s, asked, b, b / asked, NR, a[int((NR + 1) / 2)] }'
    done
  done
done
