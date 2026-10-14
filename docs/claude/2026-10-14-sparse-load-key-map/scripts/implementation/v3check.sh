#!/bin/bash
# In each x86-64-v3 clone of the lookup functions: popcnt instructions and calls of libgcc's __popcountdi2.
B=${1:-~/keysfix/new/build/protal}
objdump -d --no-show-raw-insn $B | awk '/^[0-9a-f]+ <.*>:$/ {f = $0; v3 = (f ~ /arch_x86_64_v3/) && (f ~ /LookupKmers|SparseHits|SparsePrepare|ScoreAtLevel|CountRanks/)} v3 && /popcnt/ {p[f]++} v3 && /call.*__popcountdi2/ {c[f]++} v3 {seen[f] = 1} END {for (f in seen) printf "%d popcnt, %d libgcc calls: %s\n", p[f], c[f], substr(f, 1, 90)}' | c++filt
