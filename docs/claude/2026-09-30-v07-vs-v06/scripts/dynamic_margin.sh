#!/usr/bin/env bash
# Dynamic depth-identity thresholds, tried on the benchmark's SAMs: an experimental build of the 0.7 source
# (depth_rule_patch.py: the rule from $PROTAL_DEPTH_RULE) profiles 0.7's SAMs of the deep samples again
# (--profile_only) under each rule, against the full database (v07x<tag>) and the one with species missing
# (v07hox<tag>), and of the 5M-pair sample; score.py scores them with the other runs.
set -uo pipefail
B=${BENCH:-$HOME/bench07}
SRC7=${SRC7:-$HOME/fix-build/src}
E=${EXP:-$HOME/exp-build}
T=${T:-6}
HERE=$(cd $(dirname $0) && pwd)
R=$B/runs
if [ ! -x $E/protal ]; then
  mkdir -p $E
  rsync -a --delete --exclude /build/ $SRC7/ $E/src/
  python3 $HERE/depth_rule_patch.py $E/src/src/Profiling/Profiler.h
  cmake -S $E/src -B $E/src/build -DCMAKE_BUILD_TYPE=Release > $E/cmake.log 2>&1 &&
    cmake --build $E/src/build --target protal -j $T > $E/build.log 2>&1 || { tail -20 $E/build.log; exit 1; }
  cp $E/src/build/protal $E/protal
fi
rules=("top004 top:0.04" "top008 top:0.08" "gm004 genemed:0.04" "gm006 genemed:0.06" "gm008 genemed:0.08"
       "gmz3 genemedz:3" "gmz4 genemedz:4" "gmz5 genemedz:5" "sc004k1 scaled:0.04:1" "sc004k2 scaled:0.04:2")
profile() {  # name, SAM, database
  local name=$1 sam=$2 db=$3 id=$4
  [ -f $R/$name.done ] && return
  rm -rf $R/$name; mkdir -p $R/$name
  /usr/bin/time -v -o $R/$name.time $E/protal --db $db --profile_only $sam --prefix $id -o $R/$name -t $T --no_qcmsa \
    > $R/$name.log 2>&1 && touch $R/$name.done || echo "  $name failed"
}
for entry in "${rules[@]}"; do
  tag=${entry%% *} rule=${entry#* }
  echo "$(date +%T) $tag ($rule)"
  export PROTAL_DEPTH_RULE=$rule
  for point in rl100_p500000 rl150_p500000 rl250_p500000; do
    for s in 1 2; do
      id=${point}_s_$s
      profile v07x$tag.$id $R/v07.$id/$id.sam.zst $B/tune/db07 $id
      profile v07hox$tag.$id $R/v07ho.$id/$id.sam.zst $B/tune_ho/db07 $id
    done
  done
  profile v07x$tag.deep5m_s__1 $R/v07.deep5m_s__1/deep5m_s__1.sam.zst $B/tune/db07 deep5m_s__1
done
