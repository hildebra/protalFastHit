#!/usr/bin/env bash
# The seeding profile at GTDB r226 scale, on a cluster node (experiment only; docs/claude/2026-10-12-seeding-profile).
#
#   bash seedbench_cluster.sh PROTAL_CHECKOUT DATABASE.protal R1.fq.gz R2.fq.gz OUT_DIR
#
# Needs seedbench.cpp and seedbench_build.sh beside this script, cmake and a compiler as for protal, and perf. Builds
# protal at the checkout's HEAD with -g (the same code, with debug info for perf's source lines) and seedbench beside it,
# in OUT_DIR/src and OUT_DIR/build; then, with the database:
#   1. seedbench on the first 1M pairs: 8 threads, then 32 (seedbench_8.tsv, seedbench_32.tsv): protal's seeding and
#      its variants, ns per lookup, and the blocks' best scores by size;
#   2. protal on the first 2M pairs, aligning only (--no_profile), 8 threads: perf stat with cache and TLB events
#      (perf_stat.txt), the index's huge pages once it is loaded (smaps.txt), and its seeding line (run.log);
#   3. the same run under perf record (cycles, user space): the hottest symbols (perf_report.txt) and the hottest
#      instructions of the hottest symbol with their source lines (perf_annotate.txt).
# Run it on a whole node, e.g.
#   sbatch --exclusive -c 32 --mem 100G --wrap "bash seedbench_cluster.sh ~/dev/protalFastHit DB.protal R1.fq.gz R2.fq.gz runs/seedbench"
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
checkout=$1 db=$2 r1=$3 r2=$4 out=$5
jobs=${SLURM_CPUS_PER_TASK:-$(nproc)}
mkdir -p "$out"
out=$(cd "$out" && pwd)
log=$out/seedbench_cluster.log
exec > >(tee -a "$log") 2>&1
echo "date: $(date -Iseconds); host: $(hostname); job: ${SLURM_JOB_ID:-none}; checkout: $(git -C "$checkout" rev-parse --short HEAD)"
lscpu | grep -E "Model name|^CPU\(s\)|L3 cache|NUMA node\(s\)"
cat /sys/kernel/mm/transparent_hugepage/enabled

# The build: the checkout's HEAD, Release with -g, protal and seedbench.
if [ ! -x "$out/build/seedbench" ] || [ ! -x "$out/build/protal" ]; then
  rm -rf "$out/src"; mkdir -p "$out/src"
  git -C "$checkout" archive HEAD | tar -x -C "$out/src"
  EXTRA_CMAKE="-DCMAKE_CXX_FLAGS=-g -DCMAKE_C_FLAGS=-g" bash "$here/seedbench_build.sh" "$here/seedbench.cpp" "$out/src" "$out/build" "$jobs" \
      > "$out/build.log" 2>&1 || { echo "seedbench build failed:"; tail -30 "$out/build.log"; exit 1; }
  cmake --build "$out/build" -j "$jobs" --target protal >> "$out/build.log" 2>&1 || { echo "protal build failed"; tail -30 "$out/build.log"; exit 1; }
fi
"$out/build/protal" --version 2>&1 | head -1

# The first 2M pairs, on node-local disk.
tmp=${TMPDIR:-/tmp}/seedbench_$$; mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT
zcat "$r1" | head -8000000 | gzip -1 > "$tmp/R1.fq.gz"
zcat "$r2" | head -8000000 | gzip -1 > "$tmp/R2.fq.gz"

# 1. seedbench.
for t in 8 32; do
  echo "== seedbench, 1M pairs, $t threads"
  "$out/build/seedbench" "$db" "$tmp/R1.fq.gz" "$tmp/R2.fq.gz" 1000000 "$t" 3 > "$out/seedbench_$t.tsv" 2> "$out/seedbench_$t.log"
  echo "exit $?"; cat "$out/seedbench_$t.log" "$out/seedbench_$t.tsv"
done

# 2. perf stat, the huge pages, the seeding line.
run=("$out/build/protal" --db "$db" -1 "$tmp/R1.fq.gz" -2 "$tmp/R2.fq.gz" -o "$tmp/run" -t 8 --no_profile)
echo "== perf stat"
rm -rf "$tmp/run"
perf stat -e cycles,instructions,cache-references,cache-misses,L1-dcache-loads,L1-dcache-load-misses,dTLB-loads,dTLB-load-misses \
    -o "$out/perf_stat.txt" "${run[@]}" > "$out/run.log" 2>&1 &
pid=$!
# The protal process (perf's child): its resident memory and huge pages, polled every 2 s while it runs, the reading
# with the most resident memory kept (the index loaded; protal's log is buffered, so its lines come only at the end).
best=0
while kill -0 "$pid" 2> /dev/null; do
  child=$(pgrep -P "$pid" | head -1)
  if [ -n "$child" ] && [ -r "/proc/$child/smaps_rollup" ]; then
    rss=$(awk '/^Rss:/ { print $2 }' "/proc/$child/smaps_rollup" 2> /dev/null)
    if [ -n "$rss" ] && [ "$rss" -gt "$best" ]; then
      best=$rss
      grep -E "^(Rss|AnonHugePages|ShmemPmdMapped|FilePmdMapped):" "/proc/$child/smaps_rollup" > "$out/smaps.txt" 2> /dev/null
    fi
  fi
  sleep 2
done
wait "$pid"
cat "$out/perf_stat.txt" "$out/smaps.txt"
grep -E "seeding:|Aligning reads took|Load Index took|Index in memory" "$out/run.log" | cut -c1-600
for f in "$tmp"/run/misc/*_runtime.tsv; do [ -e "$f" ] && cp "$f" "$out/run_runtime.tsv"; done

# 3. perf record: the hottest symbols, and the hottest instructions of the hottest symbol.
echo "== perf record"
rm -rf "$tmp/run"
perf record -F 1999 -e cycles:u -o "$out/perf.data" "${run[@]}" > "$out/run_record.log" 2>&1
perf report -i "$out/perf.data" --stdio --no-children --sort symbol 2> /dev/null | grep -v "^#" | grep -v "^$" | head -60 > "$out/perf_report.txt"
cat "$out/perf_report.txt"
for n in 1 2 3; do
  symbol=$(perf report -i "$out/perf.data" --stdio --no-children --sort symbol 2> /dev/null | grep -v "^#" | grep -v "^$" | sed -n "${n}p" | sed -E 's/^ *[0-9.]+% +\[\.\] //')
  [ -z "$symbol" ] && continue
  echo "== annotate: $symbol" >> "$out/perf_annotate.txt"
  # the instructions by share of the symbol's samples, the hottest 60, with the source line each comes from
  perf annotate -i "$out/perf.data" --stdio -l --no-source --symbol="$symbol" 2> /dev/null | grep -E "^ *[0-9]+\.[0-9]+" | sort -rn | head -60 >> "$out/perf_annotate.txt"
  perf annotate -i "$out/perf.data" --stdio -l --symbol="$symbol" 2> /dev/null > "$out/perf_annotate_full_$n.txt"
done
head -80 "$out/perf_annotate.txt"
echo "SEEDBENCH CLUSTER DONE"
