#!/usr/bin/env bash
# massif_summary.sh FILE.txt [SNAPSHOT] - the heap over a massif run (ms_print output: snapshot, MB) and the top
# allocation sites of one detailed snapshot (default: the last detailed one before the end).
set -u
f=$1
grep -E '^ +[0-9]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+$' "$f" | tr -d ',' |
  awk '{ printf "%s %.0f MB\n", $1, $3 / 1e6 }' | paste -sd' ' | fold -w 200
snap=${2:-}
if [ -z "$snap" ]; then
  snap=$(grep -E '^ +[0-9]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+$' "$f" | awk '{ print $1 }' | tail -2 | head -1)
fi
echo "== snapshot $snap"
awk -v s="$snap" '
  $1 == s && NF == 6 { found = 1; next }
  found && /^ +[0-9]+ +[0-9,]+ +[0-9,]+/ { exit }
  found && /->/ { print }' "$f" | grep -E '^[| ]{0,12}->' | cut -c1-230 | head -40
