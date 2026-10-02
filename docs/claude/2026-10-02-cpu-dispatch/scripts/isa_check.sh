#!/usr/bin/env bash
# Which functions of a binary hold instructions beyond baseline x86-64: AVX/AVX2 (ymm registers, v-prefixed
# SSE forms), BMI1/BMI2 (andn, bextr, blsi, blsr, blsmsk, bzhi, pdep, pext, shlx, shrx, sarx, rorx, mulx),
# LZCNT/TZCNT, FMA, F16C, MOVBE, AVX-512 (zmm, k registers). Per function: the count of such instructions.
set -uo pipefail
B=$1
objdump -d --no-show-raw-insn -M intel "$B" | awk '
  /^[0-9a-f]+ <.*>:$/ { fn = $2; next }
  {
    ins = $2
    if (ins == "") next
    beyond = 0
    if ($0 ~ /zmm|k[0-7][^a-z0-9]/ && ins ~ /^v|^k/) kind = "avx512"
    else if ($0 ~ /ymm/) kind = "avx/avx2"
    else if (ins ~ /^v/) kind = "avx(xmm)"
    else if (ins ~ /^(andn|bextr|blsi|blsr|blsmsk|bzhi|pdep|pext|shlx|shrx|sarx|rorx|mulx)$/) kind = "bmi"
    else if (ins ~ /^(lzcnt|tzcnt)$/) kind = "lzcnt/tzcnt"
    else if (ins ~ /^movbe$/) kind = "movbe"
    else next
    n[fn " " kind]++
  }
  END { for (k in n) print n[k] "\t" k }' | sort -k2,2 -k1,1nr
