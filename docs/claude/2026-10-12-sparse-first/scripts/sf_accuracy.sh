#!/usr/bin/env bash
# The calls with and without --sparse_first against the truth (docs/claude/2026-10-12-sparse-first): the test samples
# of a small real world build (BUILD/work/test/points/<point>), simulated again with the recorded commands
# (simulated.json; SIMULATE in place of the recorded binary, new folders), then profiled by each arm with BUILD's
# database and models (trained without --sparse_first), 4 threads on cores 0-3. POINTS (environment) names the points,
# default the paired-end ones (rl150_p10000 rl150_p100000 rl150_p1000000: 15 samples each); the HiFi ones are
# pb_b3000000 pb_b30000000. ARMS: "plain" or cC_tT (--sparse_first C --sparse_first_score T), default
# "plain c1_t14 c2_t14 c4_t14". protal's truth line per sample (TP, FP, FN) and records written go to WORK/calls.tsv,
# F1 by point and over all to WORK/f1.tsv, the arm's wall and CPU seconds to WORK/arms.tsv, and the failed_candidate_rate
# of true and false taxa (truth-annotated profiles) to WORK/failed_rate.tsv.
#   bash sf_accuracy.sh PROTAL SIMULATE BUILD WORK
set -uo pipefail
protal=$1 simulate=$2 build=$3 work=$4
mkdir -p "$work"
pin=(taskset -c 0-3 nice -n 5)
db=$build/protal_db/database.protal
read -r -a points <<< "${POINTS:-rl150_p10000 rl150_p100000 rl150_p1000000}"
read -r -a order <<< "${ARMS:-plain c1_t14 c2_t14 c4_t14}"
arm_args() {  # cC_tT -> --sparse_first C --sparse_first_score T
  [ "$1" = plain ] && return
  local c=${1#c}; c=${c%%_*}; local t=${1##*_t}
  echo "--sparse_first $c --sparse_first_score $t"
}

# MAP (environment): a map made otherwise (the HiFi samples of sf_hifi_sim.sh): steps 1 and 2 are left out.
if [ -n "${MAP:-}" ]; then
  points=()
  map=$MAP
fi

# 1. The test samples again: the recorded command with this simulator and new output folders.
for point in "${points[@]}"; do
  dir=$work/sim/$point
  [ -s "$dir/done" ] && continue
  rm -rf "$dir"; mkdir -p "$dir"
  mapfile -t cmd < <(sed -n '/"command": \[/,/\]/p' "$build/work/test/points/$point/simulated.json" | sed '1d;$d' |
                     sed -E 's/^ *"//; s/",?$//')
  cmd[0]=$simulate
  for i in "${!cmd[@]}"; do
    [ "${cmd[$i]}" = "-o" ] && cmd[$((i + 1))]=$dir/sim
    [ "${cmd[$i]}" = "--protal_metafile" ] && cmd[$((i + 1))]=$dir/protal
  done
  "${pin[@]}" "${cmd[@]}" > "$dir/simulate.log" 2>&1 || { echo "simulation of $point failed"; tail "$dir/simulate.log"; exit 1; }
  touch "$dir/done"
done

# 2. One map of the samples: their reads (pairs, or one file of long reads) and truth; prefix = sample.
[ -n "${MAP:-}" ] || map=$work/samples.map
[ -n "${MAP:-}" ] || {
  echo -e "#OUTPUT_DIR\t$work/arms"
  echo -e "#SAMPLEID\tPREFIX\tFIRST\tSECOND\tPROFILE_TRUTH\tREAD_TYPE"
  for point in "${points[@]}"; do
    sim=$work/sim/$point/sim
    for f in "$sim"/reads/*.fq.zst "$sim"/reads/*.fq.gz; do
      [ -e "$f" ] || continue
      case $f in *_R2.fq.*) continue ;; esac
      if [[ $f == *_R1.fq.* ]]; then
        s=$(basename "$f"); s=${s%_R1.fq.*}
        echo -e "$s\t$s\t$f\t${f/_R1.fq./_R2.fq.}\t$sim/protal_goldstd/$s.profile_truth\tpe"
      else
        s=$(basename "$f"); s=${s%.fq.*}
        echo -e "$s\t$s\t$f\t-\t$sim/protal_goldstd/$s.profile_truth\tpb"
      fi
    done
  done
} > "$map"
echo "$(($(wc -l < "$map") - 2)) samples in $map"

# 3. The arms; an arm whose log ends with protal's "Run protal took" (a finished run) is not run again (FORCE=1: it is).
for arm in "${order[@]}"; do
  if [ -z "${FORCE:-}" ] && grep -q "Run protal took" "$work/$arm.log" 2> /dev/null; then echo "$arm done before"; continue; fi
  out=$work/arms/$arm; rm -rf "$out"
  # shellcheck disable=SC2046
  "${pin[@]}" /usr/bin/time -f "TIME %e %U %S" "$protal" --db "$db" --map "$map" -o "$out" -t 4 --no_strains $(arm_args "$arm") \
      > "$work/$arm.log" 2>&1
  echo "$arm exit $?"
done
# Each arm's log: "Sample X: TP a, FP b, FN c (and ...)", "Sample X: N read pairs (or reads), ... M records written", the
# sparse-first lines (reads given up) and the run's times.
echo -e "arm\tsample\tTP\tFP\tFN\trecords\tgiven_up" > "$work/calls.tsv"
echo -e "arm\twall_s\tuser_s\tsys_s" > "$work/arms.tsv"
for arm in "${order[@]}"; do
  grep "^TIME" "$work/$arm.log" | awk -v a="$arm" '{ printf "%s\t%s\t%s\t%s\n", a, $2, $3, $4 }' >> "$work/arms.tsv"
  awk -v a="$arm" '
    /^Sample [^ ]+: [0-9]+ (read pairs|reads),/ { s = $2; sub(/:$/, "", s); n = split($0, w, " "); for (i = 1; i < n; i++) if (w[i + 1] == "records" && w[i + 2] == "written") rec[s] = w[i] }
    /^Sample [^ ]+ sparse first/ { s = $2; for (i = 3; i < NF; i++) if ($(i + 1) == "of") { g[s] = $i; break } }
    /^Sample [^ ]+: TP [0-9]+, FP [0-9]+, FN [0-9]+/ { s = $2; sub(/:$/, "", s); tp = $4; fp = $6; fn = $8; gsub(/,/, "", tp); gsub(/,/, "", fp); gsub(/[,(]/, "", fn); row[s] = tp "\t" fp "\t" fn }
    END { for (s in row) printf "%s\t%s\t%s\t%s\t%s\n", a, s, row[s], rec[s], (s in g) ? g[s] : 0 }' "$work/$arm.log" >> "$work/calls.tsv"
done

# 4. failed_candidate_rate of the true (truth 1) and false (truth 0) taxa of the truth-annotated profiles, by arm.
echo -e "arm\ttruth\ttaxa\tmean_failed_candidate_rate" > "$work/failed_rate.tsv"
for arm in "${order[@]}"; do
  find "$work/arms/$arm" -name "*.truth_annotated" -exec cat {} + 2> /dev/null | awk -F'\t' -v a="$arm" '
    $1 == "truth" { for (i = 1; i <= NF; i++) if ($i == "failed_candidate_rate") c = i; next }
    c && $1 != "truth" { n[$1]++; s[$1] += $c }
    END { for (t in n) printf "%s\t%s\t%d\t%.4f\n", a, t, n[t], s[t] / n[t] }' >> "$work/failed_rate.tsv"
done

# 5. F1 by arm and point, and over all samples.
awk -F'\t' 'NR > 1 { d = $2; sub(/_s_[0-9]+$/, "", d); tp[$1, d] += $3; fp[$1, d] += $4; fn[$1, d] += $5; rec[$1, d] += $6
              tp[$1, "all"] += $3; fp[$1, "all"] += $4; fn[$1, "all"] += $5; rec[$1, "all"] += $6; arms[$1]; ds[d]; ds["all"] }
     END { printf "arm\tpoint\tTP\tFP\tFN\tprecision\trecall\tF1\trecords\n"
           for (a in arms) for (d in ds) { p = tp[a, d] / (tp[a, d] + fp[a, d]); r = tp[a, d] / (tp[a, d] + fn[a, d])
             printf "%s\t%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\n", a, d, tp[a, d], fp[a, d], fn[a, d], p, r, 2 * p * r / (p + r), rec[a, d] } }' \
    "$work/calls.tsv" | sort -k2,2 -k1,1 | tee "$work/f1.tsv"
cat "$work/arms.tsv" "$work/failed_rate.tsv"
