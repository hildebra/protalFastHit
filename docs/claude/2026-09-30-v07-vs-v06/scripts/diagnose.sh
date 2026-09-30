#!/usr/bin/env bash
# Where 0.7's abundance bias against other strains comes from: 0.7 on the deep samples with each
# suspect turned off (runs as profile.sh's, scored by score.py like the others).
#   v07dm1.<sample>  --depth_identity_margin 1: every read counts towards abundance (0.04 by default)
#   v07wr.<sample>   --whole_read_alignment: short reads aligned as a whole, as 0.6 did
set -uo pipefail
B=${BENCH:-$HOME/bench07}
P7=${P7:-$HOME/fix-build/bin/protal}
TEST=${TEST:-$HOME/tune/V2/test/points}
T=${T:-6}
R=$B/runs
run() {
  local name=$1; shift
  [ -f $R/$name.done ] && return
  rm -rf $R/$name; mkdir -p $R/$name
  echo "$(date +%T) $name"
  /usr/bin/time -v -o $R/$name.time "$@" -o $R/$name -t $T > $R/$name.log 2>&1 && touch $R/$name.done || echo "  $name failed"
}
for point in rl100_p500000 rl150_p500000 rl250_p500000; do
  for s in 1 2; do
    id=${point}_s_$s
    r=$TEST/$point/sim/reads/$id
    pair=(-1 ${r}_R1.fq.gz -2 ${r}_R2.fq.gz --prefix $id --no_qcmsa)
    run v07dm1.$id $P7 --db $B/tune/db07 --depth_identity_margin 1 "${pair[@]}"
    run v07wr.$id $P7 --db $B/tune/db07 --whole_read_alignment "${pair[@]}"
  done
done
