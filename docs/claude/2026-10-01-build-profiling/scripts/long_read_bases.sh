#!/bin/bash
# Bases asked for and bases simulated in each long-read sample of a collection: pbsim3 makes reads for a
# sequence until their bases reach depth x its length, so every contig whose quota is a base or more gets at
# least one whole read.
# usage: long_read_bases.sh COLLECTION_DIR...
for dir in "$@"; do
  for samples in "$dir"/points/{pb,ont}_b*/sim/samples.tsv; do
    [ -f "$samples" ] || continue
    point=$(basename "$(dirname "$(dirname "$samples")")")
    asked=${point#*_b}
    tail -n +2 "$samples" | while IFS=$'\t' read -r sample reads truth community; do
      zcat "$reads" | awk -v s="$sample" -v asked="$asked" 'NR % 4 == 2 { n++; b += length($0) }
        END { printf "%-28s asked %12d bases, got %12d in %7d reads (%.1fx)\n", s, asked, b, n, b / asked }'
    done
  done
done
