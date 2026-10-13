#!/bin/bash
# Item 5: alignment threads and reading threads against the CPUs. protal (the working tree, ~/pgzfix/build/protal) on
# the mini DB (~/pgzfix/mini_db, scripts/mini_db/build_mini_db.sh), 2M read pairs (mini_reads.pl) as 64 KiB gzip members
# and as one member, on 4 CPUs (taskset -c 0-3): -t 4 with 1, 2 and 4 reading threads per file (PROTAL_INFLATE_THREADS;
# members: workers + the file's own thread, one member: in pieces), and -t 3 with 4 (a CPU left). The input line gives the
# share of the alignment spent getting reads. ROUNDS alternated rounds.
set -uo pipefail
G=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-13-gzip-members-regression/scripts
W=~/pgzfix
D=$W/threads
mkdir -p $D
cd $D
if [ ! -s members_R2.fq.gz ]; then
  zstd -dc $W/mini_db/protal_db/full_reference.fna.zst > ref.fna
  perl $G/mini_reads.pl ref.fna 2000000 5 reads
  for m in 1 2; do
    taskset -c 0-3 nice -n 5 gzip -6 -c reads_R$m.fq > one_R$m.fq.gz
    taskset -c 0-3 nice -n 5 split -b 64K --filter='gzip -6' reads_R$m.fq > members_R$m.fq.gz
    rm -f reads_R$m.fq
  done
fi
ls -la *.gz
uptime
for round in $(seq 1 ${ROUNDS:-2}); do
  for kind in members one; do
    for setting in "4 1" "4 2" "4 4" "3 4"; do
      set -- $setting
      out=out_${kind}_t$1_i$2
      rm -rf $out
      start=$(date +%s.%N)
      PROTAL_INFLATE_THREADS=$2 taskset -c 0-3 $W/build/protal --db $W/mini_db/protal_db -t $1 -1 ${kind}_R1.fq.gz -2 ${kind}_R2.fq.gz \
        --prefix s -o $out --no_strains > $out.log 2>&1
      status=$?
      wall=$(awk -v a=$start -v b=$(date +%s.%N) 'BEGIN { printf "%.2f", b - a }')
      aligning=$(grep -m1 '^Aligning reads took' $out.log | sed 's/Aligning reads took //')
      input=$(grep -m1 '^Sample s input: ' $out.log | sed 's/^Sample s input: //')
      printf 'r%s\t%s\t-t %s\tinflate %s\texit %s\twall %s\taligning %s\t%s\n' $round $kind $1 $2 $status $wall "$aligning" "$input"
    done
  done
done
uptime
