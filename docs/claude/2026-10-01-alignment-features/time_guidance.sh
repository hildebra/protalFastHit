#!/bin/bash
# Alignment time of the paired-end test samples of a build_gtdb_database.py run, alignment only (--no_profile), with:
#   A: the build without the read consensus and mate guidance (~/protal-feat)
#   B: the build with both (~/protal-cons)
#   C: the same with --no_mate_guidance
# alternated, ROUNDS rounds (the machine's wall clock varies), SAMs removed after each run; THREADS threads pinned to
# CPUS at the lowest priority, so that other work on the machine (another benchmark pinned elsewhere) is not disturbed.
#     bash time_guidance.sh RUN OUT [ROUNDS] [THREADS] [CPUS]
set -eu
RUN=$1
OUT=$2
ROUNDS=${3:-2}
THREADS=${4:-6}
CPUS=${5:-0-$(( $(nproc) - 1 ))}
mkdir -p $OUT
awk -F'\t' -v OFS='\t' '/^#/ {print; next} $8 == "pe" {print}' $RUN/test/profile_all/samples.map > $OUT/pe.map
echo "$(grep -vc '^#' $OUT/pe.map) paired-end samples"
printf "variant\tround\twall_s\tuser_s\tmax_rss_kb\n" > $OUT/times.tsv
for round in $(seq 1 $ROUNDS); do
  for variant in A B C; do
    case $variant in
      A) bin=~/protal-feat/build/protal; extra="" ;;
      B) bin=~/protal-cons/build/protal; extra="" ;;
      C) bin=~/protal-cons/build/protal; extra="--no_mate_guidance" ;;
    esac
    dir=$OUT/$variant.$round
    mkdir -p $dir/sam
    awk -F'\t' -v OFS='\t' -v O="$dir" '/^#OUTPUT_DIR/ {print "#OUTPUT_DIR", O; next} /^#/ {print; next}
        {n=split($4, p, "/"); $4 = O "/sam/" p[n]; $5 = O "/" $1; $6 = O "/" $1 ".profile"; print}' $OUT/pe.map > $dir/map
    /usr/bin/time -f "%e %U %M" -o $dir/time nice -n 19 taskset -c $CPUS $bin --db $RUN/training_db --map $dir/map \
        -t $THREADS --no_profile $extra > $dir/log 2> $dir/err
    read wall user rss < $dir/time
    printf "%s\t%s\t%s\t%s\t%s\n" $variant $round $wall $user $rss >> $OUT/times.tsv
    grep -h "fragments had one mate" $dir/log | head -1 > $dir/guidance.txt || true
    [ "$round" = "1" ] && [ "$variant" != "A" ] && cp -r $dir/sam $OUT/sam.$variant || true
    rm -rf $dir/sam
    echo "$variant round $round: ${wall}s wall, ${user}s user $(cat $dir/guidance.txt)"
  done
done
column -t $OUT/times.tsv
