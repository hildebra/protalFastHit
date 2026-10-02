#!/usr/bin/env bash
# The target_clones build as on a CPU without x86-64-v3: run under gdb, let libgcc's __cpu_indicator_init (called by
# the first ifunc resolver, before main and before any constructor) fill in __cpu_model and __cpu_features2, then
# clear their feature bits. Every resolver after that, and protal's own __builtin_cpu_supports("avx2"), sees a CPU
# without AVX2, POPCNT or SSE4 (zlib-ng, libdeflate and zstd ask the CPU themselves). Counts the calls of a few
# .default and .arch_x86_64_v3 clones, and compares the outputs with the plain x86-64 build (1,000 pairs, -t 1).
set -uo pipefail
W=$HOME/mt-work/onebcpu; rm -rf $W; mkdir -p $W
CL=${CL:-$HOME/mt-work/oneb/src/build/protal}; BASE=$HOME/mt-work/oneb-ref/src/build/protal
DB=$HOME/bench071/V071/protal_db; I=$HOME/mt-work/isacg
head -4000 $I/r1.fq > $W/r1.fq; head -4000 $I/r2.fq > $W/r2.fq
nm -S $CL | grep -E ' __cpu_model$| __cpu_features2$'
sz2=$(nm -S $CL | awk '$4 == "__cpu_features2" { print strtonum("0x" $2) }')
args="--db $DB -1 $W/r1.fq -2 $W/r2.fq --prefix s -t 1 --no_qcmsa"
# mangled clone names to count
fn() { nm $CL | awk -v p="$1" -v s="$2" '$3 ~ p && $3 ~ ("\." s "$") { print $3; exit }'; }
probes=""
for f in CigarANI GetFromLookup AlignAnchor DecodeChunk; do
  for s in default arch_x86_64_v3; do probes="$probes $(fn $f $s)"; done
done
{
  echo "set pagination off"; echo "set confirm off"
  echo "break __cpu_indicator_init"
  echo "run $args -o $W/out.nocpu > $W/log.nocpu 2>&1"
  echo "finish"
  echo "printf \"before: features %#x, features2 %#x\n\", *(unsigned int*)((char*)&__cpu_model + 12), *(unsigned int*)&__cpu_features2"
  echo "set var *(unsigned int*)((char*)&__cpu_model + 12) = 0"
  for ((i = 0; i < sz2; i += 4)); do echo "set var *(unsigned int*)((char*)&__cpu_features2 + $i) = 0"; done
  echo "delete 1"
  for p in $probes; do echo "break '$p'"; echo "commands"; echo "silent"; echo "continue"; echo "end"; done
  echo "continue"
  echo "info breakpoints"
} > $W/gdb.cmd
gdb -batch -x $W/gdb.cmd $CL > $W/gdb.log 2>&1
grep -E "before:|exited|breakpoint already hit|^[0-9]+ +breakpoint" $W/gdb.log | sed -E 's/in (protal::)?/ /' | cut -c1-200
$BASE $args -o $W/out.base > $W/log.base 2>&1
r=$(diff -r -q -x '*_runtime.tsv' -x '*.sam.zst' $W/out.base $W/out.nocpu > $W/diff && echo same || echo DIFFER)
cmp -s <(zstdcat $W/out.base/s.sam.zst) <(zstdcat $W/out.nocpu/s.sam.zst) && s="SAM text same" || s="SAM text DIFFER"
echo "without x86-64-v3 vs plain x86-64 build: outputs $r, $s"
