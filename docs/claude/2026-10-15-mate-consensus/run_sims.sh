#!/bin/bash
# The reduced r226 database's simulated paired-end test points (~/c4test/db/work/test/points) profiled on the finished
# database with one --mate_guidance mode: one protal run per point, its map holding the point's 10 samples with their
# truth (so that each profile gets a .truth_annotated table of every taxon's features and the model's call), into
# OUT/<mode>/<point>/{alignments,profiles}. The rl150 points' reads are those resim.sh made again (the build removed
# the originals after profiling them); the rl250 points' are the originals.
# DB=<database folder> profiles on another database (the build's training database: give its pe model with --model).
# LABEL=<name> names the run folder instead of the mode.
# usage: [DB=...] [POINTS="..."] [LABEL=...] run_sims.sh <protal> <mode> <out folder> [more protal options]     (WSL; 4 cores)
set -u
PROTAL=$1; MODE=$2; OUT=$3; shift 3
DB=${DB:-$HOME/c4test/db/protal_db}
P=$HOME/c4test/db/work/test/points
SIMS=${SIMS:-$HOME/mc/sims}
POINTS=${POINTS:-"rl150_p10000 rl150_p100000 rl150_p1000000 rl250_p10000 rl250_p100000"}
for point in $POINTS; do
  reads=$P/$point/sim/reads
  [ -n "$(ls $reads 2>/dev/null)" ] || reads=$SIMS/$point/sim/reads
  dir=$OUT/${LABEL:-$MODE}/$point
  # profiles/ made beforehand: protal makes it with the first profile, after that sample's truth annotation
  mkdir -p $dir/profiles $dir/alignments
  {
    printf '#OUTPUT_DIR\t%s\n#INPUT_DIR\t%s\n' "$dir" "$reads"
    printf '#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\tSAM\tPROFILE\tPROFILE_TRUTH\n'
    for t in $P/$point/sim/protal_goldstd/*.profile_truth; do
      s=$(basename $t .profile_truth)
      printf '%s\t%s\t%s_R1.fq.zst\t%s_R2.fq.zst\tpe\t%s.sam.zst\t%s.profile\t%s\n' $s $s $s $s $s $s $t
    done
  } > $dir/map.tsv
  /usr/bin/time -v taskset -c 0-3 nice -n 5 $PROTAL --db $DB --map $dir/map.tsv -t 4 --no_strains --no_qcmsa --force \
      --mate_guidance $MODE "$@" > $dir/protal.log 2> $dir/time.log
  echo "$point ${LABEL:-$MODE}: exit $? $(grep -E 'Elapsed' $dir/time.log)"
done
