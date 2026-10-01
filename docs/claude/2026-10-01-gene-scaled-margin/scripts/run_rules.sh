#!/usr/bin/env bash
# run_rules.sh WORLD_DIR - the depth rules on a world (the stress worlds of 2026-09-30, or congener_world.sh's):
#  1. the experimental build of this protal (rule_patch.py: the rule from $PROTAL_DEPTH_RULE) in $EXP;
#  2. the databases, built by this protal from the world's release, so with its gene conservation factors:
#     db7_full (every species) and db7_missing (without heldout.txt), unless present;
#  3. the samples simulated from the design's manifests and aligned once per database (--no_profile) into
#     align/<db>, unless there (the first two worlds' SAMs of 0.7 are reused: the reference is the same);
#  4. each rule (those $RULES names, by default all): --profile_only of those SAMs against db7_<db>, --knob 0
#     --no_strains, into prof7/<rule>/<db>.
# Scores: score.py WORLD_DIR.
set -uo pipefail
B=$1
SRC=${SRC:-$HOME/fix-build/src}
P=${P:-$HOME/fix-build/bin/protal}
SIM=${SIM:-$HOME/fix-build/bin/simulate_metagenomes}
E=${EXP:-$HOME/exp-build-gc}
T=${T:-6}
HERE=$(cd $(dirname $0) && pwd)
PATCH=$HERE/rule_patch.py
mkdir -p $B/logs

if [ ! -x $E/protal ] || [ $PATCH -nt $E/protal ] || [ $P -nt $E/protal ]; then
  mkdir -p $E
  rsync -a --delete --exclude /build/ $SRC/ $E/src/
  python3 $PATCH $E/src/src/Profiling/Profiler.h
  cmake -S $E/src -B $E/src/build -DCMAKE_BUILD_TYPE=Release > $E/cmake.log 2>&1 &&
    cmake --build $E/src/build --target protal -j $T > $E/build.log 2>&1 || { tail -20 $E/build.log; exit 1; }
  cp $E/src/build/protal $E/protal
fi

for db in full missing; do
  [ -f $B/db7_$db/database.protal ] && continue
  echo "$(date +%T) build db7_$db"
  extra=$([ $db = missing ] && echo "--exclude_species $B/heldout.txt" || true)
  python3 $SRC/scripts/mini_db/gtdb_to_protal_db.py --gtdb $B/gtdb --outdir $B/db7_$db -t $T $extra > $B/logs/convert7_$db.log 2>&1 &&
    $P --build --no_profile -t $T --db $B/db7_$db --reference $B/db7_$db/reference.fna \
      --full_reference $B/db7_$db/full_reference.fna --compress_level 3 > $B/logs/build7_$db.log 2>&1 || { echo "build failed"; exit 1; }
  grep "Gene conservation:" $B/logs/build7_$db.log
done

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
    echo "$(date +%T) align $s against db7_$db"
    mkdir -p $out
    r=$B/sim_rl$rl/reads/$s
    $P --db $B/db7_$db -1 ${r}_R1.fq.gz -2 ${r}_R2.fq.gz --prefix $s -o $out -t $T --no_profile > $out.log 2>&1 ||
      echo "  failed, see $out.log"
  done
done

# tag, $PROTAL_DEPTH_RULE (- for protal's own), --depth_identity_margin, --gene_conservation (- for the database's)
rules=("none - 1 -" "top004 - 0.04 none" "top008 - 0.08 none" "gtop3 - 0.08 -" "gtop0 gtop:0:0.08 0.08 -"
       "gtop2 gtop:0.02:0.08 0.08 -" "gtop4 gtop:0.04:0.08 0.08 -" "gtop3m10 - 0.10 -" "gm008 genemed:0.08 0.08 -"
       "gmedc33 gmedc:0.03:0.03 0.08 -" "gmedc32 gmedc:0.03:0.02 0.08 -"
       "gsplit13s gsplit:1.3:0.03:0.08 0.08 -" "gsplit08s gsplit:0.8:0.03:0.08 0.08 -" "gsplit07s gsplit:0.7:0.03:0.08 0.08 -"
       "gsplit08u gsplit:0.8:0.08:0.08 0.08 -"
       "gcong05a gcong:0.05:0 0.08 -" "gcong05d gcong:0.05:1 0.08 -")
for entry in "${rules[@]}"; do
  set -- $entry
  tag=$1 rule=$2 margin=$3 conservation=$4
  # $RULES: the tags to run (default all)
  [ -z "${RULES:-}" ] || [[ " $RULES " == *" $tag "* ]] || continue
  echo "$(date +%T) rule $tag"
  for db in full missing; do
    for s in $samples; do
      out=$B/prof7/$tag/$db/$s
      [ -f $out/done ] && continue
      rm -rf $out; mkdir -p $out
      args=(--depth_identity_margin $margin)
      [ $conservation = - ] || args+=(--gene_conservation $conservation)
      PROTAL_DEPTH_RULE=$([ $rule = - ] && echo "" || echo $rule) $E/protal --db $B/db7_$db \
        --profile_only $B/align/$db/$s/$s.sam.zst --prefix $s -o $out -t $T --no_qcmsa --no_strains --knob 0 \
        "${args[@]}" > $out/log 2>&1 && touch $out/done || echo "  $tag $db $s failed"
    done
  done
done
echo "$(date +%T) done"
