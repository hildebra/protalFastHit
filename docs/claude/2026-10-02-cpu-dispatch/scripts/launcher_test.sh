#!/usr/bin/env bash
# The launcher's choice with fake /proc/cpuinfo contents: a copy of protal_launcher that reads a fake file
# instead, next to stub binaries that print their name. And which of the flags this machine has.
set -uo pipefail
REPO=$(cd "$(dirname "$0")/../../../.." && pwd)
W=$(mktemp -d); trap 'rm -rf $W' EXIT
for b in protal_avx2 protal_baseline; do printf '#!/bin/sh\necho %s "$@"\n' $b > $W/$b; chmod +x $W/$b; done
sed "s#/proc/cpuinfo#$W/cpuinfo#" $REPO/protal_launcher > $W/protal; chmod +x $W/protal
v3="fpu sse sse2 ssse3 sse4_1 sse4_2 popcnt avx avx2 bmi1 bmi2 fma f16c movbe abm xsave"
try() {  # try NAME FLAGS-LINE (or "none" for no file) [ENV-ASSIGNMENT]
  rm -f $W/cpuinfo
  [ "$2" = none ] || printf 'processor\t: 0\n%s\n' "$2" > $W/cpuinfo
  printf '%-34s -> %s\n' "$1" "$(env -u PROTAL_NO_AVX2 ${3:-} $W/protal --version)"
}
try "x86-64-v3 flags" "flags		: $v3"
for missing in avx avx2 bmi1 bmi2 fma f16c movbe abm; do
  try "without $missing" "flags		: $(echo " $v3 " | sed "s/ $missing / /")"
done
try "avx2 only as part of avx512f" "flags		: fpu sse sse2 avx512f avx bmi1 bmi2 fma f16c movbe abm"
try "no /proc/cpuinfo" none
try "ARM (Features, no flags)" "Features	: fp asimd evtstrm aes pmull sha1 sha2 crc32"
try "x86-64-v3 flags, PROTAL_NO_AVX2=1" "flags		: $v3" PROTAL_NO_AVX2=1
echo "this machine's flags have: $(for f in avx avx2 bmi1 bmi2 fma f16c movbe abm; do grep -m1 '^flags' /proc/cpuinfo | grep -qw $f && printf '%s ' $f || printf '[no %s] ' $f; done)"
