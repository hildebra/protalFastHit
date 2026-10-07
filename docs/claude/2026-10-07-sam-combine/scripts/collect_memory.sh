#!/usr/bin/env bash
# collect_memory.sh - copies memory.sh's tables and the heap profiles' timelines into results/memory
# (docs/claude/2026-10-07-sam-combine).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
M=${M:-$HOME/samcombine/mem}
D=$HERE/../results/memory
mkdir -p "$D"
for tag in before after packed trim; do cp "$M/$tag.tsv" "$D/rss_$tag.tsv"; done
for m in before_strains packed_strains nostrains; do
  bash "$HERE/massif_summary.sh" "$M/massif_${m}_2.txt" 0 | head -6 > "$D/heap_${m}_2.txt"
done
bash "$HERE/evidence_size.sh" "$M/run_packed_16_yes" s1 > "$D/evidence_s1.txt"
ls -la "$D"
