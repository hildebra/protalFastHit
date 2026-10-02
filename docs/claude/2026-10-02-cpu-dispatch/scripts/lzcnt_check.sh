#!/usr/bin/env bash
# lzcnt_check.sh BINARY: lzcnt and tzcnt instructions per function, outside the AVX2/AVX-512 variants.
objdump -d --no-show-raw-insn "$1" | awk '/^[0-9a-f]+ <.*>:$/ { fn = $2; next } $2 == "lzcnt" || $2 == "tzcnt" { n[fn "\t" $2]++ }
  END { for (k in n) print n[k] "\t" k }' | c++filt | grep -v -i "avx\|vpclmul" | sed -E 's/\(.*//' | sort -k3
