#!/usr/bin/env bash
# Exact costs: callgrind instruction counts of whole runs on one thread (protal -t 1), PAR runs at a time on the
# 4 pinned cores. Alignment only (--no_profile): 0.7.8 (own database) against 0.8.1, and 0.8.1 with its new alignment
# work switched off; with profiling: 0.8.1 on the pe and PacBio subsets (--no_strains --no_qcmsa).
# Each run: callgrind.out.<name>, its log, and callgrind_annotate's inclusive top functions (<name>.top.txt).
set -uo pipefail
W=$HOME/perf9; S=$W/samples; R=${OUT:-$W/runs/callgrind}; par=${PAR:-3}
mkdir -p "$R"
pe=(-1 "$S/pe100k_R1.fq.gz" -2 "$S/pe100k_R2.fq.gz")
cat > "$R/runs.txt" << EOF
pe.v078.align|v078|$W/db078|${pe[*]} --no_profile
pe.v081.align|v081|$W/arc/protal_db|${pe[*]} --no_profile
pe.v081.align_noadapt|v081|$W/arc/protal_db|${pe[*]} --no_profile --adaptive_candidates 0
pe.v081.align_noalleles|v081|$W/arc/protal_db|${pe[*]} --no_profile --no_allele_scores
pe.v081.align_nonb|v081|$W/arc/protal_db|${pe[*]} --no_profile --no_gene_neighbours
pb.v078.align|v078|$W/db078|-1 $S/long/pb_5000000.fq.gz --read_type pb --no_profile
pb.v081.align|v081|$W/arc/protal_db|-1 $S/long/pb_5000000.fq.gz --read_type pb --no_profile
ont.v078.align|v078|$W/db078|-1 $S/long/ont_5000000.fq.gz --read_type ont --no_profile
ont.v081.align|v081|$W/arc/protal_db|-1 $S/long/ont_5000000.fq.gz --read_type ont --no_profile
pe.v081.full|v081|$W/arc/protal_db|${pe[*]} --no_strains --no_qcmsa
pb.v081.full|v081|$W/arc/protal_db|-1 $S/long/pb_5000000.fq.gz --read_type pb --no_strains --no_qcmsa
EOF
one() {
  IFS='|' read -r name bin db args <<< "$1"
  local out=$R/work_$name
  [ -s "$R/$name.top.txt" ] && return
  rm -rf "$out"; mkdir -p "$out"
  # shellcheck disable=SC2086
  taskset -c 0-3 nice -n 5 valgrind --tool=callgrind --callgrind-out-file="$R/callgrind.out.$name" \
      "$W/$bin/build/protal" --db "$db" $args -o "$out/s" -t 1 > "$R/$name.log" 2>&1
  callgrind_annotate --inclusive=yes "$R/callgrind.out.$name" 2> /dev/null | head -80 > "$R/$name.top.txt"
  rm -rf "$out"
  echo "done $name: $(grep -m1 'PROGRAM TOTALS' "$R/$name.top.txt")"
}
export -f one; export R W
grep -v '^$' "$R/runs.txt" | xargs -P "$par" -I{} bash -c 'one "$@"' _ {}
echo DONE
