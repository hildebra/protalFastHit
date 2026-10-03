#!/usr/bin/env bash
# Baseline (committed HEAD) against the packed-index build on the local worlds: the outputs must be the
# same (SAM records compared sorted: the thread order differs; every other file byte for byte, except
# the timing tables), and the peak RSS, the index's memory line and the stage times are printed.
#   ab_packed.sh [THREADS] [ROUNDS]
WORK=${WORK:-$HOME/protal-pack}
PERF=${PERF:-$HOME/protal-perf}
T=${1:-6}
ROUNDS=${2:-1}
OUT=$WORK/ab
mkdir -p "$OUT"
secs() { awk '{ t = 0; for (i = 1; i <= NF; i++) { v = $i + 0; if ($i ~ /ms$/) t += v / 1000; else if ($i ~ /m$/) t += v * 60; else if ($i ~ /s$/) t += v } printf "%.1f", t }'; }
run() {  # name db reads_dir
  local name=$1 db=$2 reads=$3 r1 r2
  [ -d "$reads/reads" ] && reads=$reads/reads  # the simulated worlds keep their reads in a subfolder
  r1=$(ls "$reads"/*_R1.fq* 2>/dev/null | head -1); r2=$(ls "$reads"/*_R2.fq* 2>/dev/null | head -1)
  [ -n "$r1" ] || { echo "$name: no reads in $reads"; return; }
  for round in $(seq 1 "$ROUNDS"); do
    for v in base new; do
      local bin=$WORK/base/build/protal; [ $v = new ] && bin=$WORK/build/protal
      rm -rf "$OUT/$name.$v"
      /usr/bin/time -f "%M" -o "$OUT/rss" "$bin" --db "$db" -1 "$r1" -2 "$r2" -o "$OUT/$name.$v" -t "$T" --no_qcmsa > "$OUT/$name.$v.stdout" 2> "$OUT/$name.$v.stderr"
      printf '%-10s %-4s load %5s s  align %6s s  run %6s s  maxrss %8d KB  %s\n' "$name" $v \
        "$(grep -m1 '^Load Index took' "$OUT/$name.$v.stdout" | sed 's/Load Index took //' | secs)" \
        "$(grep -m1 '^Aligning reads took' "$OUT/$name.$v.stdout" | sed 's/Aligning reads took //' | secs)" \
        "$(grep -m1 '^Run protal took' "$OUT/$name.$v.stdout" | sed 's/Run protal took //' | secs)" \
        "$(tail -1 "$OUT/rss")" "$(grep -m1 '^Index in memory' "$OUT/$name.$v.stdout" | cut -c1-110)"
    done
  done
  # Compare the last round's outputs.
  local same=1
  for f in $(cd "$OUT/$name.base" && find . -type f | sort); do
    local a="$OUT/$name.base/$f" b="$OUT/$name.new/$f"
    [ -f "$b" ] || { echo "  missing in new: $f"; same=0; continue; }
    case "$f" in
      *_runtime.tsv) continue ;;
      *.sam.zst) [ "$(zstd -dc "$a" | sort | md5sum)" = "$(zstd -dc "$b" | sort | md5sum)" ] || { echo "  SAM differs: $f"; same=0; } ;;
      *.sam) [ "$(sort "$a" | md5sum)" = "$(sort "$b" | md5sum)" ] || { echo "  SAM differs: $f"; same=0; } ;;
      *) cmp -s "$a" "$b" || { echo "  differs: $f"; same=0; } ;;
    esac
  done
  local n_base n_new
  n_base=$(cd "$OUT/$name.base" && find . -type f | wc -l); n_new=$(cd "$OUT/$name.new" && find . -type f | wc -l)
  [ "$n_base" = "$n_new" ] || { echo "  file counts differ: $n_base / $n_new"; same=0; }
  [ $same = 1 ] && echo "  outputs identical ($n_base files, SAM records compared sorted)"
}
for name in ${SAMPLES:-w900 mix dense_w dense_mix}; do
  db=$PERF/db900n; case $name in dense*) db=$PERF/dbdense ;; esac
  run "$name" "$db" "$PERF/reads/$name"
done
