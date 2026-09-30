#!/usr/bin/env bash
# Species the database lacks: the tuning world's database without the species the training database
# of ~/tune/V2 leaves out (whole clades and 20% of the species), built by both versions, and the deep
# and a shallow test sample profiled against it. Reads of the species left out land on relatives, which
# the abundance of 0.7 is meant to resist (--depth_identity_margin 0.04) and 1 lets through.
#   v06ho.<sample>, v07ho.<sample>, v07hodm1.<sample> (--depth_identity_margin 1)
set -uo pipefail
B=${BENCH:-$HOME/bench07}
SRC7=${SRC7:-$HOME/fix-build/src}
P6=${P6:-$HOME/protal-0.6.0a/src/build/protal} P7=${P7:-$HOME/fix-build/bin/protal}
TEST=${TEST:-$HOME/tune/V2/test/points}
HELDOUT=${HELDOUT:-$HOME/tune/V2/heldout_species.txt}
T=${T:-6}
W=$B/tune_ho R=$B/runs
if [ ! -f $W/db07/database.protal ]; then
  rm -rf $W; mkdir -p $W
  python3 $SRC7/scripts/mini_db/gtdb_to_protal_db.py --from_db $B/tune/conv --exclude_species $HELDOUT --outdir $W/conv
  cp -r $W/conv $W/db07; cp -r $W/conv $W/db06; mv $W/db06/model_pe.xml $W/db06/model.xml
  for v in 07 06; do
    p=$([ $v = 07 ] && echo $P7 || echo $P6)
    $p --build --no_profile -t $T --db $W/db$v --reference $W/db$v/reference.fna --full_reference $W/db$v/full_reference.fna > $B/logs/tune_ho.build$v.log 2>&1
  done
  $P7 --add_model $HOME/tune/V2/trained_model.xml --read_type pe --db $W/db07 -t $T > $B/logs/tune_ho.add_model.log 2>&1
fi
run() {
  local name=$1; shift
  [ -f $R/$name.done ] && return
  rm -rf $R/$name; mkdir -p $R/$name
  echo "$(date +%T) $name"
  /usr/bin/time -v -o $R/$name.time "$@" -o $R/$name -t $T > $R/$name.log 2>&1 && touch $R/$name.done || echo "  $name failed"
}
for point in rl100_p500000 rl150_p500000 rl250_p500000 rl150_p10000; do
  for s in 1 2; do
    id=${point}_s_$s
    r=$TEST/$point/sim/reads/$id
    pair=(-1 ${r}_R1.fq.gz -2 ${r}_R2.fq.gz --prefix $id --no_qcmsa)
    run v06ho.$id $P6 --db $W/db06 "${pair[@]}"
    run v07ho.$id $P7 --db $W/db07 "${pair[@]}"
    run v07hodm1.$id $P7 --db $W/db07 --depth_identity_margin 1 "${pair[@]}"
  done
done
