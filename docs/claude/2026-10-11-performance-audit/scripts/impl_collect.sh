#!/usr/bin/env bash
# The follow-up's results and scripts into the report folder (results/impl, scripts/impl_*).
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
D=/mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-11-performance-audit
W=$HOME/impl/perf
O=$D/results/impl; mkdir -p "$O"
cp "$W/runs2/runs.tsv" "$O/runs2.tsv"; cp "$W/runs2/cohort.tsv" "$O/cohort2.tsv"
cp "$W/runs3/runs.tsv" "$O/runs3.tsv"; cp "$W/runs3/cohort.tsv" "$O/cohort3.tsv"
cp "$W/tables/profile_only.tsv" "$O/profile_only_before_header_fix.tsv"
cp "$W/tables3/profile_only.tsv" "$O/profile_only.tsv"
cp "$W/tables/callgrind.tsv" "$O/callgrind_loaders_before_header_fix.tsv"
cp "$W/post/callgrind.tsv" "$O/callgrind_loaders.tsv"
cp "$W/cg_align/callgrind.tsv" "$O/callgrind_alignment.tsv"
grep -v "^\[" "$A/impl_identity.out" > "$O/identity.txt"
grep -E "^(same|different) SAM" "$A/impl_post.out" > "$O/sam_checks.txt"
sed -n "/] identity/,/MERGED DONE/p" "$A/impl_merged.out" | grep -E "^(same|DIFFERENT)" > "$O/identity_merged_0dea1ff.txt"
grep -E "^(Wrote|Rewrite)" "$W/db_bin.log" > "$O/compress_db.txt"
for f in impl_build.sh impl_unit.sh impl_verify_src.sh impl_identity.sh impl_e2e.sh impl_perf_setup.sh impl_perf_runs.sh impl_after.sh \
         impl_tables.sh impl_cg_align.sh impl_final.sh impl_round2.sh impl_post.sh impl_asan.sh impl_runs3.sh impl_merged.sh impl_collect.sh; do
  cp "$A/$f" "$D/scripts/$f"
done
ls -la "$O" "$D/scripts" | grep -E "impl|tsv|txt"
