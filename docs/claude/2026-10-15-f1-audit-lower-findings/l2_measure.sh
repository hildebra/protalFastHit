#!/bin/bash
# F1 audit L2 (--gene_end_clips): the runs of the measurement, one after the other (4 cores each), all with
# --mate_guidance consensus (the reduced r226 database records no run settings, so auto would be guide):
#  mocks_base:   the six paired-end mocks (../2026-10-15-mate-consensus/run_mocks.sh) with the base binary (68316e7);
#  mocks_off, mocks_credit: the same with this branch's binary, --gene_end_clips off and credit;
#  fin_off, fin_credit: the simulated test points (run_sims.sh; the rl150 reads remade by resim.sh into SIMS) on the
#    finished reduced r226 database;
#  tdb_off, tdb_credit: the same on the build's training database (its held-out species absent) with the finished
#    database's pe model, as the build scored its test set;
#  mocks_taxon, fin_taxon, tdb_taxon (not by default): credit with a binary built with l2_taxon_common.patch (every
#    candidate of a taxon credited as the taxon's least clipped one), given as the fourth argument.
# usage: [STEPS="..."] [SIMS=...] l2_measure.sh <protal> <base protal> <out folder> [patched protal]   (WSL; 4 cores)
set -u
BIN=$1; BASE=$2; R=$3; TAXON=${4:-}
HERE=$(cd "$(dirname "$0")" && pwd)
D=$HERE/../2026-10-15-mate-consensus
export SIMS=${SIMS:-$HOME/f1l2/sims}
STEPS=${STEPS:-"mocks_base mocks_off mocks_credit fin_off fin_credit tdb_off tdb_credit"}
for step in $STEPS; do
  echo "== $step $(date +%T)"
  case $step in
    mocks_base)   LABEL=base bash $D/run_mocks.sh $BASE consensus $R/mocks ;;
    mocks_off)    LABEL=off bash $D/run_mocks.sh $BIN consensus $R/mocks --gene_end_clips off ;;
    mocks_credit) LABEL=credit bash $D/run_mocks.sh $BIN consensus $R/mocks --gene_end_clips credit ;;
    fin_off)      LABEL=off bash $D/run_sims.sh $BIN consensus $R/fin --gene_end_clips off ;;
    fin_credit)   LABEL=credit bash $D/run_sims.sh $BIN consensus $R/fin --gene_end_clips credit ;;
    tdb_off)      DB=$HOME/c4test/db/work/training_db LABEL=off bash $D/run_sims.sh $BIN consensus $R/tdb \
                      --model $HOME/c4test/db/model_logs/trained_model.xml --gene_end_clips off ;;
    tdb_credit)   DB=$HOME/c4test/db/work/training_db LABEL=credit bash $D/run_sims.sh $BIN consensus $R/tdb \
                      --model $HOME/c4test/db/model_logs/trained_model.xml --gene_end_clips credit ;;
    mocks_taxon)  PROTAL_L2_TAXON_COMMON=1 LABEL=taxon bash $D/run_mocks.sh $TAXON consensus $R/mocks --gene_end_clips credit ;;
    fin_taxon)    PROTAL_L2_TAXON_COMMON=1 LABEL=taxon bash $D/run_sims.sh $TAXON consensus $R/fin --gene_end_clips credit ;;
    tdb_taxon)    PROTAL_L2_TAXON_COMMON=1 DB=$HOME/c4test/db/work/training_db LABEL=taxon bash $D/run_sims.sh $TAXON consensus \
                      $R/tdb --model $HOME/c4test/db/model_logs/trained_model.xml --gene_end_clips credit ;;
  esac
done
echo "== measure done $(date +%T)"
