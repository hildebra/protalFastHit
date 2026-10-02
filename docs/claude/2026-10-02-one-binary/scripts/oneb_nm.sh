#!/usr/bin/env bash
# The clones in the target_clones build: per marked function, its .default, .arch_x86_64_v3 and .resolver
# symbols, and whether any .default clone (or anything outside the v3 clones and the known AVX2 kernels) has
# ymm instructions.
B=$HOME/mt-work/oneb/src/build/protal
nm -C $B | grep -E '\.(default|arch_x86_64_v3|resolver)' | sed -E 's/^[0-9a-f]+ . //' | sed -E 's/\(.*\)(\.[a-z_0-9]+)$/\1/' | sort | uniq -c | sort -k2 | head -60
echo "clone symbols: $(nm $B | grep -c '\.arch_x86_64_v3$') v3, $(nm $B | grep -c '\.default$') default, $(nm $B | grep -c '\.resolver$') resolvers"
# functions with ymm instructions
objdump -d --no-show-raw-insn $B | awk '/^[0-9a-f]+ <.*>:$/ { f = $2 } /ymm/ { n[f]++ } END { for (f in n) print n[f], f }' | c++filt | sort -k2 > /tmp/ymm_funcs.txt
echo "functions with ymm: $(wc -l < /tmp/ymm_funcs.txt); of them .default clones: $(grep -c '\.default>' /tmp/ymm_funcs.txt)"
grep -v 'arch_x86_64_v3>' /tmp/ymm_funcs.txt | grep -iv 'avx2\|avx\|zng\|deflate\|ZSTD\|HUF_\|FSE_\|adler\|crc32\|chunk\|compare256\|longest_match\|slide_hash\|inflate_fast\|Bases16\|FillAvx2\|ScanWindowsAvx2\|PackAvx2\|Unpack' | head -30
