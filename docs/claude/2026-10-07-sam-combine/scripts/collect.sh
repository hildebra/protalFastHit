#!/usr/bin/env bash
# collect.sh - copies run.sh's and no_outdir.sh's summaries and maps into results/ (docs/claude/2026-10-07-sam-combine).
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
W=${W:-$HOME/samcombine/work}
D=$HERE/../results
mkdir -p "$D"
cp "$W/summary.txt" "$D/run_summary.txt"
cp "$W/merged.map" "$W/sams.map" "$W/sams2.map" "$W/grow.map" "$D/"
for name in profile_only_one reads_one reads_two reads_two_no_strains reads_two_dot; do
  f=$W/no_o/$name.log
  abort=$(grep -E 'what\(\)' "$f" | head -1 | sed 's/^ *//')
  echo "$name: ${abort:-no abort}"
done > "$D/no_outdir.txt"
ls -la "$D"
cat "$D/no_outdir.txt"
