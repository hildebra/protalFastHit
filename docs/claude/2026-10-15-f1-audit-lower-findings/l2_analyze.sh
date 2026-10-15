#!/bin/bash
# F1 audit L2 (--gene_end_clips): the tables of the report from l2_measure.sh's runs, into RESULTS:
#  identical_mocks.txt: the base binary against off, per mock sample: the records (sorted, without the header), the
#    header (without the settings and @PG lines), the profile (without its # lines) and the .truth_annotated table;
#  records_<set>.tsv: off against credit by ../2026-10-15-mate-consensus/compare_records.py (fin and tdb with the
#    simulation's truth per fragment; the mocks without); changed_<set>.tsv/.txt: the simulated fragments whose
#    counted species changed, by kind (changed_fragments.py); mates_<set>.tsv: the mates whose MAPQ or taxon changed,
#    and those that crossed MAPQ 4, each way (l2_mapq_cross.py);
#  calls_<set>.txt: the calls, credit against off, by ../2026-10-15-mock-near-misses/spike_calls.py.
# B=<label> compares off with that run instead of credit (e.g. taxon, l2_taxon_common.patch); its tables are named
# <table>_<set>_<label>.
# usage: [STEPS="identical mocks fin tdb"] [B=credit] l2_analyze.sh <runs folder> <results folder>     (WSL)
set -u
R=$1; RES=$2
B=${B:-credit}
SUF=""; [ $B != credit ] && SUF=_$B
HERE=$(cd "$(dirname "$0")" && pwd)
D=$HERE/..
PY=$HOME/micromamba/envs/protal-db-build/bin/python
FIN=$HOME/c4test/db/protal_db
TDB=$HOME/c4test/db/work/training_db
WORK=$HOME/c4test/db/work
mkdir -p $RES
same() { [ "$1" = "$2" ] && echo same || echo DIFFER; }
STEPS=${STEPS:-"identical mocks fin tdb"}
for step in $STEPS; do
  case $step in
    identical)
      {
        for sam in $R/mocks/base/alignments/*.sam.zst; do
          s=$(basename $sam .sam.zst)
          off=$R/mocks/off/alignments/$s.sam.zst
          n=$(zstdcat $sam | grep -vc '^@')
          a=$(zstdcat $sam | grep -v '^@' | LC_ALL=C sort | md5sum | cut -c1-32)
          b=$(zstdcat $off | grep -v '^@' | LC_ALL=C sort | md5sum | cut -c1-32)
          ha=$(zstdcat $sam | grep '^@' | grep -v -e '^@CO.protal settings: ' -e '^@PG' | md5sum | cut -c1-32)
          hb=$(zstdcat $off | grep '^@' | grep -v -e '^@CO.protal settings: ' -e '^@PG' | md5sum | cut -c1-32)
          pa=$(grep -v '^#' $R/mocks/base/profiles/$s.profile | md5sum | cut -c1-32)
          pb=$(grep -v '^#' $R/mocks/off/profiles/$s.profile | md5sum | cut -c1-32)
          ta=$(md5sum < $R/mocks/base/profiles/$s.profile.truth_annotated | cut -c1-32)
          tb=$(md5sum < $R/mocks/off/profiles/$s.profile.truth_annotated | cut -c1-32)
          echo -e "$s\trecords $n $(same $a $b)\theader $(same $ha $hb)\tprofile $(same $pa $pb)\ttruth_annotated $(same $ta $tb)"
        done
        echo "settings lines (base, off):"
        zstdcat $(ls $R/mocks/base/alignments/*.sam.zst | head -1) | grep -m1 '^@CO.protal settings'
        zstdcat $(ls $R/mocks/off/alignments/*.sam.zst | head -1) | grep -m1 '^@CO.protal settings'
      } > $RES/identical_mocks.txt
      cat $RES/identical_mocks.txt ;;
    mocks)
      specs=()
      for sam in $R/mocks/off/alignments/*.sam.zst; do
        s=$(basename $sam .sam.zst)
        specs+=("$s:$sam:$R/mocks/$B/alignments/$s.sam.zst")
      done
      $PY -I $D/2026-10-15-mate-consensus/compare_records.py $RES/records_mocks$SUF.tsv off $B $FIN "${specs[@]}" 2> /dev/null
      $PY -I $HERE/l2_mapq_cross.py $RES/mates_mocks$SUF.tsv "${specs[@]}" 2> /dev/null
      $PY -I $D/2026-10-15-mock-near-misses/spike_calls.py $R/mocks/$B/profiles $R/mocks/off/profiles > $RES/calls_mocks$SUF.txt ;;
    fin|tdb)
      dbdir=$FIN; [ $step = tdb ] && dbdir=$TDB
      specs=()
      mkdir -p $RES/tables_$step/off $RES/tables_$step/$B
      for dir in $R/$step/off/*/; do
        point=$(basename $dir)
        for sam in $dir/alignments/*.sam.zst; do
          s=$(basename $sam .sam.zst)
          specs+=("$s:$sam:$R/$step/$B/$point/alignments/$s.sam.zst")
          ln -sf $dir/profiles/$s.profile.truth_annotated $RES/tables_$step/off/
          ln -sf $R/$step/$B/$point/profiles/$s.profile.truth_annotated $RES/tables_$step/$B/
        done
      done
      $PY -I $D/2026-10-15-mate-consensus/compare_records.py $RES/records_$step$SUF.tsv off $B $dbdir --work $WORK "${specs[@]}" 2> /dev/null
      $PY -I $D/2026-10-15-mate-consensus/changed_fragments.py $RES/changed_$step$SUF.tsv $dbdir $WORK "${specs[@]}" > $RES/changed_$step$SUF.txt
      $PY -I $HERE/l2_mapq_cross.py $RES/mates_$step$SUF.tsv "${specs[@]}" 2> /dev/null
      $PY -I $D/2026-10-15-mock-near-misses/spike_calls.py $RES/tables_$step/$B $RES/tables_$step/off > $RES/calls_$step$SUF.txt
      rm -rf $RES/tables_$step ;;
  esac
done
echo "analysis done"
