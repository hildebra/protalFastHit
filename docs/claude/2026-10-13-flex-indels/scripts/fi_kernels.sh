#!/usr/bin/env bash
# The flex scan's instructions per step without and with --flex_indels' scores (docs/claude/2026-10-13-flex-indels), from
# the compiled kernels in the unit-test binary (where the tests call them, not inlined): each loop of BestAvx2<false/true>
# (8 cells a step) and ScanAvx512<false/true> (16 cells a step; its second loop is the ties'), as its instructions from
# the jump back to its start, found by objdump. The AVX-512 counts are of the code, which runs only on an AVX-512 CPU.
#   bash fi_kernels.sh PROTAL_TESTS
set -uo pipefail
binary=$1
objdump -d -C --no-show-raw-insn "$binary" | awk '
  /^[0-9a-f]+ <.*>:$/ { name = $0; keep = (name ~ /flex_scan::(BestAvx2|ScanAvx512|ScoreScalar)<(true|false)>/); n = 0; next }
  keep && /^ *[0-9a-f]+:\t/ {
    addr = $1; sub(/:$/, "", addr); a = strtonum("0x" addr); line[++n] = a; text[n] = $0
    if ($0 ~ /\tj[a-z]+ +[0-9a-f]+ </) {
      target = $0; sub(/.*\tj[a-z]+ +/, "", target); sub(/ <.*/, "", target); t = strtonum("0x" target)
      if (t < a) {  # a jump back: the loop from t to here
        count = 0; simd = 0
        for (i = 1; i <= n; i++) if (line[i] >= t) { count++; if (text[i] ~ /vpopcnt|vpshufb|vpmaddubsw|vpternlog|vplzcnt/) simd++ }
        fn = name; sub(/^[0-9a-f]+ </, "", fn); sub(/\(.*$/, "", fn)
        if (count > 6) printf "%s\tloop of %d instructions (%d popcount/shuffle/ternary-logic)\n", fn, count, simd
      }
    }
  }'
