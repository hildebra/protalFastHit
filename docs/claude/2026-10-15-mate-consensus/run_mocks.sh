#!/bin/bash
# The six paired-end mock samples in ~/Mocks (D6331 DNA, D6331 HKU and Auckland, Treichel, hCom2, hCom1a) profiled on
# the reduced r226 database with one --mate_guidance mode, in one protal run whose map gives each sample its truth's
# internal taxids (PROFILE_TRUTH: the profile gets a .truth_annotated table of every taxon's features), into
# OUT/<mode>/{alignments,profiles}; the profiles are also linked into OUT/score/mc_<mode>/<set>/profiles/, the layout
# scripts/mock_benchmark/score_mocks.py scores.
# LABEL=<name> names the run folder instead of the mode.
# usage: [LABEL=...] run_mocks.sh <protal> <mode> <out folder> [more protal options]     (WSL; 4 cores)
set -u
PROTAL=$1; MODE=$2; OUT=$3; shift 3
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$HERE/../../..
DB=$HOME/c4test/db/protal_db
T=$REPO/scripts/mock_benchmark/truth
PY=$HOME/micromamba/envs/protal-db-build/bin/python
READS=$HOME/Mocks/reads
dir=$OUT/${LABEL:-$MODE}
mkdir -p $dir/truth $dir/profiles $dir/alignments  # profiles/ before the first sample's truth annotation
# set, sample, truth table (scripts/mock_benchmark/samples.tsv and the local test's samples_extra.tsv)
SAMPLES="d6331_dna:d6331_dna_SRR38158090:d6331 d6331_pe_labs:d6331_hku_SRR23926922:d6331
d6331_pe_labs:d6331_auckland_SRR33710519:d6331 treichel:treichel_ERR14788576:treichel_even70
hcom2_invitro:hcom2_SRR21819229:hcom2_closed hcom1a_dropouts:hcom1a_SRR24799833:hcom1a_full"
{
  printf '#OUTPUT_DIR\t%s\n' "$dir"
  printf '#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\tSAM\tPROFILE\tPROFILE_TRUTH\n'
  for spec in $SAMPLES; do
    IFS=: read -r set sample truth <<< "$spec"
    $PY -I $HERE/../2026-10-13-ecoli-hifi/truth_taxids.py $T/$truth.tsv $DB/genome2tiid.tsv > $dir/truth/$sample.truth 2> /dev/null
    first=$(awk -F'\t' -v s=$sample '$2==s{print $4}' $READS/downloaded.tsv)
    second=$(awk -F'\t' -v s=$sample '$2==s{print $5}' $READS/downloaded.tsv)
    printf '%s\t%s\t%s\t%s\tpe\t%s.sam.zst\t%s.profile\t%s\n' $sample $sample $first $second $sample $sample $dir/truth/$sample.truth
  done
} > $dir/map.tsv
/usr/bin/time -v taskset -c 0-3 nice -n 5 $PROTAL --db $DB --map $dir/map.tsv -t 4 --no_strains --no_qcmsa --force \
    --mate_guidance $MODE "$@" > $dir/protal.log 2> $dir/time.log
echo "mocks ${LABEL:-$MODE}: exit $? $(grep -E 'Elapsed' $dir/time.log)"
for spec in $SAMPLES; do
  IFS=: read -r set sample truth <<< "$spec"
  mkdir -p $OUT/score/mc_${LABEL:-$MODE}/$set/profiles
  ln -sf $dir/profiles/$sample.profile $OUT/score/mc_${LABEL:-$MODE}/$set/profiles/$sample.profile
done
