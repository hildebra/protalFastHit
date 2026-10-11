#!/usr/bin/env bash
# The audit's measurements in order, after the archaea world's finished database is built (build_arcdb.sh): the models
# (the small real world's pe and PacBio models, trained by 0.8.x on 170 species, stored for all four read types),
# 0.7.8's own database, then m3 (start-up tables), m1 (alignment A/B), m5 (profiling), m4 (cohort memory) and
# m2 (callgrind). Every step pinned to cores 0-3.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
W=$HOME/perf9; P=$W/v081/build/protal; DB=$W/arc/protal_db; M=$HOME/realworld_s1/model_logs
step() { echo "[$(date +%T)] $*"; }
until grep -q "^exit" "$A/build_arcdb.out" 2> /dev/null; do sleep 20; done
step "finished database: $(grep '^exit' "$A/build_arcdb.out")"
for t in pe:trained_model.xml se:trained_model.xml pb:trained_model_pb.xml ont:trained_model_pb.xml; do
    "$P" --add_model "$M/${t#*:}" --read_type "${t%%:*}" --db "$DB" > "$W/add_model_${t%%:*}.log" 2>&1
    step "model ${t%%:*}: exit $? $(tail -1 "$W/add_model_${t%%:*}.log" | cut -c1-150)"
done
step "db078"; bash "$A/db078.sh" > "$W/db078.out" 2>&1; step "$(grep -E 'exit|FAIL' "$W/db078.out" | head -2)"
step "m3 tables"; bash "$A/m3_tables.sh" > "$W/m3.out" 2>&1; step "$(tail -1 "$W/m3.out")"
step "m1 align"; bash "$A/m1_align.sh" > "$W/m1.out" 2>&1; step "$(tail -1 "$W/m1.out")"
step "m5 profiling"; bash "$A/m5_profiling.sh" > "$W/m5.out" 2>&1; step "$(tail -1 "$W/m5.out")"
step "m4 cohort"; bash "$A/m4_cohort.sh" > "$W/m4.out" 2>&1; step "$(tail -1 "$W/m4.out")"
step "m2 callgrind"; PAR=4 bash "$A/m2_callgrind.sh" > "$W/m2.out" 2>&1; step "$(tail -1 "$W/m2.out")"
step "ALL DONE"
