#!/bin/bash
# The speed tests of items 3, 4 and 5 on a quiet machine, at a758a4a (git archive -o a758a4a.tar a758a4a, from Git
# Bash, next to this script; not kept). Each timed run starts only when no other protal or python3 runs and cores 0-3
# were under 10% busy for 5 s (/proc/stat); all runs on cores 0-3, in alternated rounds; medians at the end.
#   A (item 4): bench_pgz on one gzip member in pieces, the look-ahead before a758a4a (V0: gap 2T + 4, ahead gap + 2T)
#     against a758a4a's (V3: gap T + 2, ahead gap + T), at 2, 3 and 4 threads: wall seconds, max RSS.
#   B (item 3): bench_pgz (V3) with the read-ahead hints and with PROTAL_READ_AHEAD=0, the file dropped from the page
#     cache first (cold) or not (warm): one member in pieces on 4 threads and on one, 64 KiB members on 4.
#   C (item 5): protal on the mini DB, 2M read pairs as one member and as 64 KiB members, -t 4 with 1, 2 and 4 reading
#     threads per file (PROTAL_INFLATE_THREADS): aligning seconds and the reader seconds a thread.
set -uo pipefail
HERE=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-reading-input
G=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-gzip-members-regression/scripts
W=~/pgzfix
B=~/pgzread
OUT=$HERE/results
mkdir -p $W/src $B/data $W/threads
quiet() {  # no other protal or python3, and cores 0-3 under 10% busy over 5 s
  local a b
  while true; do
    if ! pgrep -x protal > /dev/null && ! pgrep -x python3 > /dev/null; then
      a=$(awk '/^cpu[0-3] / { u += $2 + $3 + $4 + $7 + $8; t += $2 + $3 + $4 + $5 + $6 + $7 + $8 } END { print u, t }' /proc/stat)
      sleep 5
      b=$(awk '/^cpu[0-3] / { u += $2 + $3 + $4 + $7 + $8; t += $2 + $3 + $4 + $5 + $6 + $7 + $8 } END { print u, t }' /proc/stat)
      awk -v a="$a" -v b="$b" 'BEGIN { split(a, x, " "); split(b, y, " "); exit !((y[1] - x[1]) < 0.1 * (y[2] - x[2])) }' && return 0
    fi
    sleep 10
  done
}
# Medians of column $3 (and $4 if given) of the rows in $1 by the key columns $2 (e.g. "2,3").
medians() {
  awk -F'\t' -v keys="$2" -v col="$3" -v col2="${4:-0}" '
    BEGIN { nk = split(keys, kc, ",") }
    { k = ""; for (i = 1; i <= nk; i++) k = k (i > 1 ? "\t" : "") $kc[i]; n = ++c[k]; v[k, n] = $col; if (col2) w[k, n] = $col2 }
    function med(arr, key, n,   i, j, t, x) { for (i = 1; i <= n; i++) x[i] = arr[key, i] + 0
        for (i = 1; i <= n; i++) for (j = i + 1; j <= n; j++) if (x[j] < x[i]) { t = x[i]; x[i] = x[j]; x[j] = t }
        lo = x[1]; hi = x[n]; return n % 2 ? x[(n + 1) / 2] : (x[n / 2] + x[n / 2 + 1]) / 2 }
    END { for (k in c) { m = med(v, k, c[k]); l = lo; h = hi; printf "%s\t%.3f\t(%.3f-%.3f, %d runs)", k, m, l, h, c[k]
                         if (col2) printf "\t%.2f", med(w, k, c[k]); printf "\n" } }' "$1" | sort
}

# --- build a758a4a: protal, simulate_metagenomes, bench_pgz V0 and V3 ---------------------------------------------
tar -xf $HERE/scripts/a758a4a.tar -C $W/src || exit 1
taskset -c 0-3 nice -n 5 cmake -S $W/src -B $W/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=OFF \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache > $W/configure.log 2>&1 || { tail -20 $W/configure.log; exit 1; }
taskset -c 0-3 nice -n 5 ninja -C $W/build -j4 protal simulate_metagenomes > $W/build.log 2>&1 || { grep error: $W/build.log | head; exit 1; }
S=$W/src
cp $S/src/IO/ParallelGzip.cpp $S/src/IO/ParallelGzip_V3.cpp
sed -e 's/if (impl.gap == 0) impl.gap = size_t{impl.threads} + 2;/if (impl.gap == 0) impl.gap = 2 * size_t{impl.threads} + 4;/' \
    -e 's/impl.ahead = impl.gap + size_t{impl.threads};/impl.ahead = impl.gap + 2 * size_t{impl.threads};/' \
    $S/src/IO/ParallelGzip.cpp > $S/src/IO/ParallelGzip_V0.cpp
grep -q 'impl.gap = 2 \* size_t{impl.threads} + 4;' $S/src/IO/ParallelGzip_V0.cpp && grep -q 'impl.ahead = impl.gap + 2 \* size_t' $S/src/IO/ParallelGzip_V0.cpp || { echo "V0 not patched"; exit 1; }
for v in V0 V3; do
  taskset -c 0-3 nice -n 5 g++ -O3 -std=c++20 -march=x86-64 -mtune=generic -pthread -DNDEBUG -DPGZ_NEW -I $S/src -I $S/src/Utilities \
    -I $S/lib $G/bench_pgz.cpp $S/src/IO/ParallelGzip_$v.cpp -lisal -lzstd -o $B/bench_$v || exit 1
done
echo "built: $($W/build/protal --version 2>&1 | head -1)"

# --- data ------------------------------------------------------------------------------------------------------------
cd $B/data
if [ ! -s m64k.gz ]; then
  taskset -c 0-3 nice -n 5 perl $G/make_fastq.pl 4000000 7 > R1.fq
  taskset -c 0-3 nice -n 5 gzip -6 -c R1.fq > one.gz
  taskset -c 0-3 nice -n 5 split -b 64K --filter='gzip -6' R1.fq > m64k.gz
  rm -f R1.fq
fi
if [ ! -f $W/mini_db/protal_db/database.protal ]; then
  (cd $S && PROTAL=$W/build/protal taskset -c 0-3 nice -n 5 bash scripts/mini_db/build_mini_db.sh $W/mini_db > $W/mini_db.log 2>&1) || { echo "mini db failed"; exit 1; }
fi
cd $W/threads
if [ ! -s members_R2.fq.gz ]; then
  zstd -dc $W/mini_db/protal_db/full_reference.fna.zst > ref.fna
  perl $G/mini_reads.pl ref.fna 2000000 5 reads
  for m in 1 2; do
    taskset -c 0-3 nice -n 5 gzip -6 -c reads_R$m.fq > one_R$m.fq.gz
    taskset -c 0-3 nice -n 5 split -b 64K --filter='gzip -6' reads_R$m.fq > members_R$m.fq.gz
    rm -f reads_R$m.fq
  done
fi
echo "data ready"; uptime

# --- A: the look-ahead ---------------------------------------------------------------------------------------------
cd $B/data
: > $OUT/speed_lookahead.tsv
for round in $(seq 1 ${ROUNDS:-5}); do
  for t in 2 3 4; do
    for v in V0 V3; do
      quiet
      line=$(taskset -c 0-3 $B/bench_$v one.gz $t)
      printf 'r%s\t%s\t%s\t%s\t%s\n' $round $v $t "$(grep -oE 'wall [0-9.]+' <<< "$line" | cut -d' ' -f2)" \
        "$(grep -oE 'maxrss_mb [0-9]+' <<< "$line" | cut -d' ' -f2)" >> $OUT/speed_lookahead.tsv
    done
  done
done
echo "A done"

# --- B: read-ahead, cold and warm ----------------------------------------------------------------------------------
: > $OUT/speed_readahead.tsv
for round in $(seq 1 ${ROUNDS:-5}); do
  for case in "one.gz 4" "one.gz 1" "m64k.gz 4"; do
    set -- $case
    for hints in 1 0; do
      for cache in cold warm; do
        quiet
        [ $cache = cold ] && dd if=$1 iflag=nocache count=0 status=none
        line=$(PROTAL_READ_AHEAD=$hints taskset -c 0-3 $B/bench_V3 $1 $2)
        printf 'r%s\t%s\t%s\t%s\t%s\t%s\n' $round $1 $2 $hints $cache "$(grep -oE 'wall [0-9.]+' <<< "$line" | cut -d' ' -f2)" \
          >> $OUT/speed_readahead.tsv
      done
    done
  done
done
echo "B done"

# --- C: reading threads at -t 4 on 4 CPUs --------------------------------------------------------------------------
cd $W/threads
: > $OUT/speed_threads.tsv
seconds() { awk '{ s = 0; for (k = 1; k <= NF; k++) { v = $k + 0; if ($k ~ /ms$/) s += v / 1000; else if ($k ~ /m$/) s += 60 * v; else if ($k ~ /s$/) s += v } printf "%.3f", s }' <<< "$1"; }
for round in $(seq 1 ${ROUNDS_C:-6}); do
  for kind in one members; do
    for i in 1 2 4; do
      quiet
      rm -rf out_c
      PROTAL_INFLATE_THREADS=$i taskset -c 0-3 $W/build/protal --db $W/mini_db/protal_db -t 4 -1 ${kind}_R1.fq.gz -2 ${kind}_R2.fq.gz \
        --prefix s -o out_c --no_strains > out_c.log 2>&1
      a=$(seconds "$(grep -m1 '^Aligning reads took' out_c.log | sed 's/Aligning reads took //')")
      s=$(grep -m1 '^Sample s input: ' out_c.log | grep -oE 'spent [0-9.]+ s each' | grep -oE '[0-9.]+')
      printf 'r%s\t%s\t%s\t%s\t%s\n' $round $kind $i "$a" "${s:-NA}" >> $OUT/speed_threads.tsv
    done
  done
done
grep -m1 '^Sample s input: ' out_c.log > $OUT/speed_threads_input_line.txt
echo "C done"; uptime

echo; echo "A: look-ahead (variant, threads): median wall s (range, runs), median max RSS MB"
medians $OUT/speed_lookahead.tsv 2,3 4 5
echo; echo "B: read-ahead (file, threads, hints, cache): median wall s"
medians $OUT/speed_readahead.tsv 2,3,4,5 6
echo; echo "C: threads (input, reading threads per file): median aligning s, median reader s a thread"
medians $OUT/speed_threads.tsv 2,3 4 5
