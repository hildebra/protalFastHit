#!/usr/bin/env bash
# --depth_identity_margin swept on the deep samples, by profiling 0.7's SAMs again (--profile_only),
# against the full database (v07dmNNN) and the one without the held-out species (v07hodmNNN);
# NNN = the margin in hundredths (004 = 0.04, the default).
set -uo pipefail
B=${BENCH:-$HOME/bench07}
P7=${P7:-$HOME/fix-build/bin/protal}
T=${T:-6}
R=$B/runs
for m in 0.04 0.06 0.08 0.10 0.15 0.25; do
  tag=$(printf "%03d" $(echo "$m * 100 / 1" | bc))
  for point in rl100_p500000 rl150_p500000 rl250_p500000; do
    for s in 1 2; do
      id=${point}_s_$s
      for db in full ho; do
        if [ $db = full ]; then name=v07dm$tag.$id src=$R/v07.$id dbdir=$B/tune/db07
        else name=v07hodm$tag.$id src=$R/v07ho.$id dbdir=$B/tune_ho/db07; fi
        [ -f $R/$name.done ] && continue
        rm -rf $R/$name; mkdir -p $R/$name
        sam=$(ls $src/$id.sam* | head -1)
        /usr/bin/time -v -o $R/$name.time $P7 --db $dbdir --profile_only $sam --prefix $id -o $R/$name -t $T \
          --depth_identity_margin $m --no_qcmsa > $R/$name.log 2>&1 && touch $R/$name.done || echo "  $name failed"
      done
    done
  done
done
