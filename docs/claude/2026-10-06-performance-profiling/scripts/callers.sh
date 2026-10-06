#!/usr/bin/env bash
# callers.sh <callgrind.out> <regex>...: the callers of each matching function (callgrind_annotate --tree=caller).
f=$1; shift
callgrind_annotate --tree=caller --inclusive=yes "$f" 2>/dev/null > /tmp/perf6_callers.txt
for re in "$@"; do
  echo "=== $re"
  # A caller block: lines starting with spaces and '<' then the function line ('*' marks it).
  awk -v re="$re" '
    /^ *[0-9,]+ .*< / { buf = buf $0 "\n"; next }
    /\*  / { if ($0 ~ re) { printf "%s%s\n\n", buf, $0 } ; buf = ""; next }
    { buf = "" }
  ' /tmp/perf6_callers.txt | cut -c1-230 | head -${LINES_MAX:-40}
done
