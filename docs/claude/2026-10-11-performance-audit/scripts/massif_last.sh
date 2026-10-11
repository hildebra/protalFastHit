#!/usr/bin/env bash
# The heap's totals per snapshot, and the allocation tree of the largest late snapshot (two levels; names shortened).
M=$HOME/perf9/runs/massif/massif.txt
awk '/^ *[0-9]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+$/ { print $1, $3 }' "$M" | tail -30 | tr '\n' ';'; echo
last=$(awk '/^ *[0-9]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+$/ { n = $1 } END { print n }' "$M")
pick=${SNAP:-$((last - 1))}
echo "== snapshot $pick"
awk -v s="$pick" '
  /^ *[0-9]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+ +[0-9,]+$/ { cur = $1 }
  cur == s && (/^->/ || /^\| ->/ || /^\|   ->/) { print }' "$M" | sed -E 's/\(in [^)]*\)//; s/std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >/string/g' | cut -c1-230 | head -70
