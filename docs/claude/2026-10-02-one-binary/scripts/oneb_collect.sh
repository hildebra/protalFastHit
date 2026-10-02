#!/usr/bin/env bash
# Results of the one-binary work into the report folder: callgrind totals, per-function differences, clone
# symbols, the gdb runs of both paths, the ASan run.
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/1680d2f3-dfad-4cd7-9113-3d5aefceb3b5/scratchpad
R=/mnt/c/Users/hildebra/Documents/locDev/protal-clones/docs/claude/2026-10-02-one-binary/results
mkdir -p $R
W=$HOME/mt-work/onebcg
{ printf "workload\tbuild\tinstructions\tvs_base_pct\n"
  for w in pe100k ont3M prof500k; do
    b0=$(grep 'PROGRAM TOTALS' $W/$w.base.txt | awk '{ gsub(",", "", $1); print $1 }')
    for b in base v2 v3 cl cl2 cl3; do
      n=$(grep 'PROGRAM TOTALS' $W/$w.$b.txt | awk '{ gsub(",", "", $1); print $1 }')
      printf "%s\t%s\t%s\t%.2f\n" $w $b $n $(echo "100 * ($n - $b0) / $b0" | bc -l)
    done
  done; } > $R/callgrind_totals.tsv
for w in pe100k ont3M prof500k; do
  { echo "# $w: self instructions per function (clone suffixes merged), the 25 largest savings"
    for p in "base v2" "base v3" "base cl2" "cl2 base" "v3 cl2" "cl2 v3"; do
      set -- $p; echo; echo "## $1 -> $2 (positive: $2 executes fewer)"; A=$1 B=$2 bash $S/oneb_diff.sh $w 25 2>&1 | cut -c1-220
    done; } > $R/functions_$w.txt
done
B=$HOME/mt-work/oneb/src/build/protal
{ echo "# clone symbols in the target_clones build (protal), demangled"; nm -C $B | grep -E '\[clone \.(default|arch_x86_64_v3|resolver)\]$' | sed -E 's/^[0-9a-f]+ . //' | sort
  echo; echo "clone symbols: $(nm $B | grep -c '\.arch_x86_64_v3$') v3, $(nm $B | grep -c '\.default$') default, $(nm $B | grep -c '\.resolver$') resolvers"
  objdump -d --no-show-raw-insn $B | awk '/^[0-9a-f]+ <.*>:$/ { f = $2 } /ymm/ { n[f]++ } END { for (f in n) print n[f], f }' | c++filt | sort -k2 > /tmp/ymm.txt
  echo "functions with ymm registers: $(wc -l < /tmp/ymm.txt), of them .default clones: $(grep -c '\.default\]>' /tmp/ymm.txt)"
  echo; echo "# functions with ymm registers (count, name)"; sed -E 's/\(.*\)//' /tmp/ymm.txt; } > $R/clones.txt
C=$HOME/mt-work/onebcpu
{ echo "# feature bits cleared after libgcc's first CPU check (oneb_nocpu.sh): breakpoint hits"; grep -E "before:|exited normally|already hit|^[0-9]+ +breakpoint" $C/gdb.log | sed -E 's/\(.*\[clone/ [clone/' | cut -c1-160
  echo; echo "# feature bits as the CPU reports them (oneb_v3cpu.sh): breakpoint hits"; grep -E "before:|exited normally|already hit|^[0-9]+ +breakpoint" $C/gdb_v3.log | sed -E 's/\(.*\[clone/ [clone/' | cut -c1-160; } > $R/dispatch_gdb.txt
{ cat $S/oneb_asan2.out; grep -E "tests passed|Total Test time" $HOME/mt-work/oneb-asan/ctest.log; } > $R/asan.txt
{ cat $S/oneb_cg.out $S/oneb_rebuild.out $S/oneb_final.out; grep "vs base" $S/oneb_time.out | sed 's/^/timed runs, round 1: /'; } | grep -E "vs base" > $R/outputs_identical.txt
grep -vE '^\[Thread|breakpoint +keep' $S/oneb_final2.out > $R/final_checks.txt
ls -la $R
