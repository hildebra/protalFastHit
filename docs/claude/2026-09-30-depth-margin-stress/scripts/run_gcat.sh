#!/usr/bin/env bash
# The second stress world, whose genes differ in conservation (simulate_gtdb_release.py --gene_rates
# categories), and the rules that scale per gene: world.sh and run.sh (the rules of the first world) into
# $STRESS (default ~/stress2); then each gene's conservation factor, estimated from the database's genes as
# protal could (estimate_gene_rates.py: per gene, and per category) and the simulator's true rates; then the
# per-gene rules with each of them ($PROTAL_GENE_RATES). Scores: score_stress.py $STRESS.
set -uo pipefail
export STRESS=${STRESS:-$HOME/stress2}
B=$STRESS
SRC7=${SRC7:-$HOME/fix-build/src}
E=${EXP:-$HOME/exp-build}
T=${T:-6}
HERE=$(cd $(dirname $0) && pwd)
RELEASE_ARGS="--gene_rates categories" bash $HERE/world.sh || exit 1
bash $HERE/run.sh || exit 1

# Conservation factors: the release converted again (the build packed reference.fna into database.protal).
if [ ! -f $B/rates/estimated.gene.tsv ]; then
  mkdir -p $B/rates
  python3 $SRC7/scripts/mini_db/gtdb_to_protal_db.py --gtdb $B/gtdb --outdir $B/conv_full -t $T > $B/logs/convert_rates.log 2>&1
  python3 $HERE/estimate_gene_rates.py $B/conv_full $SRC7/scripts/mini_db/markers_r226.tsv $B/rates/estimated \
    $B/gtdb/simulation/gene_rates.tsv | tee $B/rates/estimate.log
fi
samples=$(cut -f1 $B/design/design.tsv | tail -n +2 | sort -u)
for source in true gene category; do
  rates=$B/rates/estimated.$source.tsv
  for entry in "gtop3k5 gtop:0.03:0.05" "gtop4k4 gtop:0.04:0.04" "gmedc3k2 gmedc:0.03:0.02" "gmedc3k3 gmedc:0.03:0.03"; do
    tag=${entry%% *}_$source rule=${entry#* }
    echo "$(date +%T) rule $tag ($rule, $rates)"
    for db in full missing; do
      for s in $samples; do
        out=$B/prof/$tag/$db/$s
        [ -f $out/done ] && continue
        rm -rf $out; mkdir -p $out
        PROTAL_DEPTH_RULE=$rule PROTAL_GENE_RATES=$rates $E/protal --db $B/db_$db \
          --profile_only $B/align/$db/$s/$s.sam.zst --prefix $s -o $out -t $T --no_qcmsa --no_strains --knob 0 \
          --depth_identity_margin 0.08 > $out/log 2>&1 && touch $out/done || echo "  $tag $db $s failed"
      done
    done
  done
done
echo "$(date +%T) done"
