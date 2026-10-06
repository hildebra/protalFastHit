#!/usr/bin/env bash
# annot.sh <callgrind.out> <file suffix> [n]: the n costliest lines of that source file (inlined code included),
# from callgrind_annotate's auto-annotation of the -g build. The full annotation is cached next to the .out file.
f=$1; suffix=$2; n=${3:-30}
A=${f%.out}.annot.txt
[ -s "$A" ] || callgrind_annotate --auto=yes --context=0 "$f" > "$A" 2>/dev/null
awk -v suf="$suffix" '
  /^-- Auto-annotated source:/ { on = (index($0, suf) > 0); next }
  on && !/=>/ && /^ *[0-9][0-9,]* \(/ { line = $0; c = $1; gsub(",", "", c); print c "\t" line }
' "$A" | sort -t$'\t' -k1,1 -n -r | head -$n | cut -f2- | cut -c1-190
