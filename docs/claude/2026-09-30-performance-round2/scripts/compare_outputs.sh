#!/bin/bash
# Byte-identical outputs of two builds on the same reads: compare_outputs.sh BIN_A BIN_B DB READS_DIR THREADS [extra args]
# Runs the whole pipeline (alignment, profile; qcmsa as asked) with each binary into its own folder and compares every
# output file except the ones that hold timings (*_runtime.tsv). With THREADS > 1 the SAM records come in a
# different order from run to run, so SAMs are compared sorted (and, with one thread, as files too).
# SE=1 aligns the first file only, as single-end reads.
source "$(dirname "$0")/env.sh"
A=$1; B=$2; db=$3; reads=$4; t=$5; shift 5
d=$PERF_DIR/runs2/compare; rm -rf $d; mkdir -p $d/a $d/b
r1=$(ls $reads/*_R1.fq* | head -1); r2=$(ls $reads/*_R2.fq* | head -1)
second=(-2 $r2); [ -n "${SE:-}" ] && second=()
for x in a b; do
  bin=$A; [ $x = b ] && bin=$B
  $bin --db $db -1 $r1 "${second[@]}" -o $d/$x/out -t $t --no_qcmsa "$@" > $d/$x.stdout 2> $d/$x.stderr
  rc=$?
  echo "$x: rc=$rc $(grep -cE 'alignments' $d/$x.stdout) lines about alignments in stdout"
  [ $rc = 0 ] || { echo "protal failed: $(head -c 300 $d/$x.stderr)"; exit 1; }
done
fail=0
n=0
while IFS= read -r f; do
  rel=${f#$d/a/}
  case $rel in *_runtime.tsv) continue ;; esac
  n=$((n+1))
  if [ ! -e $d/b/$rel ]; then echo "MISSING in b: $rel"; fail=1; continue; fi
  if cmp -s $f $d/b/$rel; then continue; fi
  case $rel in
    *.sam.zst|*.sam.gz|*.sam)
      dec() { case $1 in *.zst) zstd -dc $1 ;; *.gz) zcat $1 ;; *) cat $1 ;; esac; }
      if [ "$(dec $f | sort | md5sum)" = "$(dec $d/b/$rel | sort | md5sum)" ]; then
        [ $t = 1 ] && { echo "SAM bytes differ (records equal when sorted): $rel"; fail=1; } || echo "sam equal sorted: $rel"
      else echo "SAM RECORDS DIFFER: $rel"; fail=1; fi ;;
    *) echo "DIFFERENT: $rel"; fail=1 ;;
  esac
done < <(find $d/a -type f | sort)
[ $n -gt 0 ] || fail=1
echo "compared $n files (without *_runtime.tsv): $([ $fail = 0 ] && echo IDENTICAL || echo DIFFERENCES)"
