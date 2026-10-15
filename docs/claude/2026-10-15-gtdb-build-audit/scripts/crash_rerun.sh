#!/bin/bash
# usage: crash_rerun.sh CHECKOUT
# r226 v28 in small: a build that stops at the parity check after protal_db's build in the background has ended, then a
# rerun into the same OUTDIR with the fixed scripts and a new, empty scratch folder (a new cluster job's node-local
# disk). The first run uses the parity check before b42cbc3, made to wait 90 s so that protal_db is finished first (in
# v28 it was finished 1 h 45 min before the check), with setups whose first and last design point by name have no
# single-end samples (rl100_p*, rl75_p*), and a few simulations streamed.
HERE=$(cd "$(dirname "$0")" && pwd)
CHECKOUT=$1
rm -rf ~/bgaudit/src_v28 ~/bgaudit/probe/crash2 ~/bgaudit/probe/crash2_scratch ~/bgaudit/probe/crash2_scratch2
cp -a ~/bgaudit/src ~/bgaudit/src_v28
git -C "$CHECKOUT" show b42cbc3^:scripts/check_model_parity.py |
  sed 's/^    opts = parse_args(argv)$/    opts = parse_args(argv)\n    import time; time.sleep(90)/' \
  > ~/bgaudit/src_v28/scripts/check_model_parity.py
OPTS=(--read-setups 75:HS25:250:40,100:HS20:300:40 --se-read-setups 150:HSXt:350:50 --stream-above 0.0005)
PROBE_SRC=~/bgaudit/src_v28/scripts/mini_db bash "$HERE/run_probe.sh" crash2 "${OPTS[@]}" | head -1
cp ~/bgaudit/probe/crash2/console.log ~/bgaudit/probe/crash2_console_run1.log
grep -E "Command failed" ~/bgaudit/probe_crash2.out | cut -c1-200
bash "$HERE/run_probe.sh" crash2 "${OPTS[@]}" --scratch ~/bgaudit/probe/crash2_scratch2 | head -1
echo "== the rerun: what it kept and what it made again"
awk '/^--- /{n++} n==2' ~/bgaudit/probe/crash2/console.log | grep -E "protal_db|training_db|kept|again|earlier run" |
  grep -v "so far" | cut -c1-200
echo "== marker truth: lines on the console (first run, rerun), the rerun's rows with marker reads"
grep -c "Marker truth" ~/bgaudit/probe/crash2_console_run1.log
awk '/^--- /{n++} n==2' ~/bgaudit/probe/crash2/console.log | grep -c "Marker truth"
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="meta_marker_reads") c=i; next} {n++; if($c!="" && $c!="-") f++}
  END{print "rows", n, "with marker reads", f+0}' ~/bgaudit/probe/crash2_scratch2/training/training_data.tsv
