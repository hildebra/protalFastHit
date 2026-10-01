#!/usr/bin/env bash
# Speed of parallel profiling: the 0.7.1 binary (base) against the build of this change (par), alternated, niced.
#  A: profiling only (--profile_only) of the 5M-pair SAM at -t 1, 2, 3, 4, 6, three times each
#  B: profiling only of the eight full-database SAMs together at -t 6, twice each
#  C: whole runs (alignment and profiling) of the 5M-pair sample at -t 6, twice each
# One row per run in ~/mt-work/speed/runs.tsv: name, binary, threads, rep, profiling s, wall s, peak RSS kB, load.
set -uo pipefail
W=$HOME/mt-work/speed
B=$HOME/bench071
DB=$B/V071/protal_db
SAM=$HOME/mt-audit/runs/b6/s1.sam.zst
R1=$B/samples_deep/points/rl150_p5000000/sim/reads/rl150_p5000000_s_1_R1.fq.gz
R2=${R1%_R1.fq.gz}_R2.fq.gz
mkdir -p $W
[ -f $W/runs.tsv ] || printf "set\tbinary\tthreads\trep\tprofiling_s\talign_s\twall_s\tpeak_rss_kb\tload\n" > $W/runs.tsv
seconds() {  # "1m 2s 345ms" / "2s 3ms" / "345ms" -> seconds
  awk '{ s = 0; for (i = 1; i <= NF; i++) { v = $i; if (v ~ /ms$/) { sub(/ms/, "", v); s += v / 1000 } else if (v ~ /m$/) { sub(/m/, "", v); s += v * 60 } else if (v ~ /s$/) { sub(/s/, "", v); s += v } } printf "%.3f", s }'
}
run() {  # run SET BIN THREADS REP ARGS...
  local set=$1 bin=$2 t=$3 rep=$4; shift 4
  local d=$W/$set.$bin.t$t.r$rep
  rm -rf $d; mkdir -p $d
  local load=$(cut -d' ' -f1 /proc/loadavg)
  nice -n 5 /usr/bin/time -v -o $d.time $HOME/mt-work/bin/protal-$bin --db $DB "$@" -o $d/out -t $t --no_qcmsa > $d.log 2>&1
  local prof=$(grep 'Profiling took' $d.log | sed 's/.*took //' | seconds)
  local align=$(grep 'Aligning reads took' $d.log | sed 's/.*took //' | seconds)
  local wall=$(grep 'Elapsed' $d.time | awk '{print $NF}' | awk -F: '{ if (NF == 3) print $1 * 3600 + $2 * 60 + $3; else print $1 * 60 + $2 }')
  local rss=$(grep 'Maximum resident' $d.time | awk '{print $NF}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $set $bin $t $rep "$prof" "${align:-}" "$wall" "$rss" "$load" | tee -a $W/runs.tsv
  rm -rf $d/out/*.sam.zst  # keep the profiles, not the SAMs
}
ALL=""; PREFIXES=""
for s in rl150_p500000 rl100_p500000 rl150_p10000 rl150_p1000; do
  ALL+="${ALL:+,}$B/runs/v071.full.pe.${s}_s_1/${s}_s_1.sam.zst"; PREFIXES+="${PREFIXES:+,}pe_$s"
done
ALL+=",$B/runs/v071.full.se.rl150_p500000_s_1/rl150_p500000_s_1.sam.zst,$B/runs/v071.full.pb.pb_b90000000_s_1/pb_b90000000_s_1.sam.zst"
ALL+=",$B/runs/v071.full.ont.ont_b90000000_s_1/ont_b90000000_s_1.sam.zst,$SAM"
PREFIXES+=",se,pb,ont,pe5M"
for rep in 1 2 3; do
  for bin in base par2; do run B2 $bin 6 $rep --profile_only $ALL --prefix $PREFIXES; done
done
echo done
