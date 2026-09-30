#!/usr/bin/env bash
# The stress test's runs (after world.sh): the samples simulated from the design's manifests; aligned once
# per database by 0.7 (--no_profile); then profiled again from those SAMs (--profile_only, --knob 0 so that
# every taxon with reads is reported) under each depth-identity rule, by the experimental build of the
# benchmark (../../2026-09-30-v07-vs-v06/scripts/depth_rule_patch.py; the rule from $PROTAL_DEPTH_RULE).
set -uo pipefail
B=${STRESS:-$HOME/stress}
SRC7=${SRC7:-$HOME/fix-build/src}
P7=${P7:-$HOME/fix-build/bin/protal}
SIM=${SIM:-$HOME/fix-build/bin/simulate_metagenomes}
E=${EXP:-$HOME/exp-build}
T=${T:-6}
HERE=$(cd $(dirname $0) && pwd)
PATCH=$HERE/../../2026-09-30-v07-vs-v06/scripts/depth_rule_patch.py

# The experimental build, again when the patch is newer than it.
if [ ! -x $E/protal ] || [ $PATCH -nt $E/protal ]; then
  mkdir -p $E
  rsync -a --delete --exclude /build/ $SRC7/ $E/src/
  python3 $PATCH $E/src/src/Profiling/Profiler.h
  cmake -S $E/src -B $E/src/build -DCMAKE_BUILD_TYPE=Release > $E/cmake.log 2>&1 &&
    cmake --build $E/src/build --target protal -j $T > $E/build.log 2>&1 || { tail -20 $E/build.log; exit 1; }
  cp $E/src/build/protal $E/protal
fi

for setup in "100 HS20 300 40" "150 HS25 350 50"; do
  set -- $setup
  rl=$1
  if [ ! -f $B/sim_rl$rl/protal.meta ]; then
    echo "$(date +%T) simulate 2x$rl"
    $SIM --from_manifest $B/design/rl$rl.manifest.tsv -o $B/sim_rl$rl --read_length $rl --sequencer $2 \
      --fragment_mean $3 --fragment_stdev $4 -t $T --protal_metafile $B/prot_rl$rl > $B/logs/simulate_rl$rl.log 2>&1 ||
      { echo "simulation failed"; exit 1; }
  fi
done

samples=$(cut -f1 $B/design/design.tsv | tail -n +2 | sort -u)
for db in full missing; do
  for s in $samples; do
    rl=${s%%_*}; rl=${rl#s}
    out=$B/align/$db/$s
    [ -f $out/$s.sam.zst ] && continue
    echo "$(date +%T) align $s against db_$db"
    mkdir -p $out
    r=$B/sim_rl$rl/reads/$s
    $P7 --db $B/db_$db -1 ${r}_R1.fq.gz -2 ${r}_R2.fq.gz --prefix $s -o $out -t $T --no_profile > $out.log 2>&1 ||
      echo "  failed, see $out.log"
  done
done

rules=("none -" "top004 top:0.04" "top008 top:0.08" "top012 top:0.12" "gm004 genemed:0.04" "gm006 genemed:0.06"
       "gm008 genemed:0.08" "gq20m004 geneq:0.2:0.04" "gq20m006 geneq:0.2:0.06" "fl006f008 floor:0.06:0.08"
       "sc004k2 scaled:0.04:2")
for entry in "${rules[@]}"; do
  tag=${entry%% *} rule=${entry#* }
  echo "$(date +%T) rule $tag ($rule)"
  for db in full missing; do
    for s in $samples; do
      out=$B/prof/$tag/$db/$s
      [ -f $out/done ] && continue
      rm -rf $out; mkdir -p $out
      margin=$([ $tag = none ] && echo 1 || echo 0.04)
      PROTAL_DEPTH_RULE=$([ $tag = none ] && echo "" || echo $rule) $E/protal --db $B/db_$db \
        --profile_only $B/align/$db/$s/$s.sam.zst --prefix $s -o $out -t $T --no_qcmsa --no_strains --knob 0 \
        --depth_identity_margin $margin > $out/log 2>&1 && touch $out/done || echo "  $tag $db $s failed"
    done
  done
done
echo "$(date +%T) done"
