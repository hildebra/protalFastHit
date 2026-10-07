#!/usr/bin/env bash
# check_fixes.sh - the fixes of docs/claude/2026-10-07-sam-combine against run.sh's work folder: runs without -o,
# --profile_only patterns (quoted and unquoted, a folder per sample), protal_map_utils merge reusing SAMs, and sample
# IDs that cannot name files. PROTAL: the fixed build; SRC: its source tree (for protal_map_utils).
set -u
PROTAL=${PROTAL:-$HOME/samcombine/src/build/protal}
SRC=${SRC:-$HOME/samcombine/src}
DB=${DB:-$HOME/samcombine/db/database.protal}
W=${W:-$HOME/samcombine/work}
RUN="taskset -c 0-3 nice -n 5"
MAPU="python3 $SRC/scripts/protal_map_utils"
cd "$W" || exit 1
rm -rf fix; mkdir -p fix/logs
SUMMARY=$W/fix/summary.txt; : > "$SUMMARY"
note() { echo "$*" | tee -a "$SUMMARY"; }
p() {
  local name=$1; shift
  (cd "${CWD:-$W}" && $RUN "$PROTAL" "$@") > "$W/fix/logs/$name.log" 2>&1
  note "$name: exit $?"
}
show() { grep -a -E "${2:-Error|Note|Warning|what\(\)|All alignments|Load index}" "$W/fix/logs/$1.log" | head -${3:-8} | sed 's/^/    | /' | tee -a "$SUMMARY"; }
rows() {
  for f in "$1"/*.raw.msa.fna; do
    [ -e "$f" ] || { note "    (no raw MSAs in $1)"; return; }
    note "    $(basename "$f" .raw.msa.fna): $(grep '>' "$f" | tr -d '>' | tr '\n' ' ')"
  done
}
same() {
  local diffs=0 n=0
  for f in "$1"/*.raw.msa.fna "$1"/*.raw.partition.txt "$1"/*.meta.tsv; do
    [ -e "$f" ] || continue
    n=$((n + 1)); cmp -s "$f" "$2/$(basename "$f")" || { diffs=$((diffs + 1)); note "    differs: $(basename "$f")"; }
  done
  note "    $n strain files of $1 against $2: $diffs differ"
}
# Same MSAs up to the row names: sequences in row order
same_seqs() {
  local diffs=0
  for f in "$1"/*.raw.msa.fna; do
    cmp -s <(grep -v '>' "$f") <(grep -v '>' "$2/$(basename "$f")") || diffs=$((diffs + 1))
  done
  note "    sequences of $1 against $2: $diffs MSAs differ"
}

note "== 1. no -o: -1/-2, two samples, in an empty folder"
mkdir -p fix/no_o_reads
CWD=$W/fix/no_o_reads p no_o_reads --db "$DB" -t 4 --no_qcmsa -1 "$W/reads/run1/sa_1.fq,$W/reads/run1/sb_1.fq" \
  -2 "$W/reads/run1/sa_2.fq,$W/reads/run1/sb_2.fq" --prefix sa,sb
note "    files: $(cd fix/no_o_reads && ls | tr '\n' ' ')"
same run1/strains fix/no_o_reads/strains

note "== 2. no -o: --profile_only of copies of four SAMs"
mkdir -p fix/no_o_po/sams fix/no_o_po/cwd
for s in run1/alignments/sa run1/alignments/sb run2/alignments/sc run2/alignments/sd; do cp "$s.sam.zst" fix/no_o_po/sams/; done
CWD=$W/fix/no_o_po/cwd p no_o_profile_only --db "$DB" -t 4 --no_qcmsa --profile_only "$W/fix/no_o_po/sams/sa.sam.zst,$W/fix/no_o_po/sams/sb.sam.zst,$W/fix/no_o_po/sams/sc.sam.zst,$W/fix/no_o_po/sams/sd.sam.zst"
note "    cwd: $(ls fix/no_o_po/cwd | tr '\n' ' '); sams: $(ls fix/no_o_po/sams | grep -c profile$) profiles next to the SAMs"
same joint/strains fix/no_o_po/cwd/strains

note "== 3. quoted pattern over run1 and run2 (distinct names)"
p pattern_quoted --db "$DB" -t 4 --no_qcmsa --profile_only 'run[12]/alignments/*.sam.zst' -o fix/pattern
show pattern_quoted "Note"
rows fix/pattern/strains
same joint/strains fix/pattern/strains

note "== 4. unquoted: the shell's expansion as arguments"
p pattern_unquoted --db "$DB" -t 4 --no_qcmsa -o fix/unquoted --profile_only run1/alignments/*.sam.zst run2/alignments/*.sam.zst
rows fix/unquoted/strains
same joint/strains fix/unquoted/strains

note "== 5. a pattern over three runs, sa in two of them: named by their folders"
p pattern_collide --db "$DB" -t 4 --no_qcmsa --profile_only 'run*/alignments/*.sam.zst' -o fix/collide
show pattern_collide "Note"
rows fix/collide/strains
note "    profiles: $(ls fix/collide | grep '\.profile$' | tr '\n' ' ')"

note "== 6. a folder per sample, every SAM named aln.sam.zst"
for s in sa sb sc sd; do
  r=run1; [ "$s" = sc ] || [ "$s" = sd ] && r=run2
  mkdir -p "fix/persample/$s"; cp "$r/alignments/$s.sam.zst" "fix/persample/$s/aln.sam.zst"
done
p pattern_per_sample --db "$DB" -t 4 --no_qcmsa --profile_only 'fix/persample/*/aln.sam.zst' -o fix/per_sample
show pattern_per_sample "Note"
rows fix/per_sample/strains
same joint/strains fix/per_sample/strains

note "== 7. errors: a pattern without SAMs, a folder, a missing file, all at once"
p pattern_errors --db "$DB" -t 4 --no_qcmsa --profile_only 'nothing/*.sam.zst,run1/alignments,run1/alignments/*.err,missing.sam.zst' -o fix/errors
show pattern_errors "Error"
p no_sam_at_all --db "$DB" -t 4 --no_qcmsa --profile_only '' -o fix/none
show no_sam_at_all "Error"

note "== 8. a stray argument in an alignment run"
p stray_argument --db "$DB" -t 4 --no_qcmsa -1 reads/run1/sa_1.fq reads/run1/sb_1.fq -2 reads/run1/sa_2.fq --prefix sa -o fix/stray
show stray_argument "rror|An option"

note "== 9. protal_map_utils merge of run1.map and run2.map: the runs' SAMs are kept"
$MAPU merge --map run1.map run2.map --out "$W/fix/merged" > fix/merged.map; note "    merge: exit $?"
sed 's/^/    | /' fix/merged.map | tee -a "$SUMMARY"
$MAPU validate --map fix/merged.map; note "    validate: exit $?"
mv reads reads.away
$MAPU validate --map fix/merged.map; note "    validate without the reads: exit $?"
p merged_map --map fix/merged.map --db "$DB" -t 4 --no_qcmsa
show merged_map "All alignments|Load index|Error"
same joint/strains fix/merged/strains
mv reads.away reads
$MAPU merge --map run1.map run2.map --out "$W/fix/merged_new" --new-sams > fix/merged_new.map; note "    merge --new-sams: exit $?"
grep -v '^#' fix/merged_new.map | cut -f1,4 | sed 's/^/    | /' | tee -a "$SUMMARY"

note "== 10. merge of single-end maps (SECOND '-'; and no SECOND column)"
mkdir -p fix/se
printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\nx\t%s\t-\tx\n' "$W/fix/se/a" "$W/reads/run1/sa_1.fq" > fix/se/a.map
printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\ny\t%s\t-\ty\n' "$W/fix/se/b" "$W/reads/run1/sb_1.fq" > fix/se/b.map
$MAPU merge --map fix/se/a.map fix/se/b.map --out "$W/fix/se/m" > fix/se/merged.map; note "    merge: exit $?"
sed 's/^/    | /' fix/se/merged.map | tee -a "$SUMMARY"
printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tPREFIX\nx\t%s\tx\n' "$W/fix/se/c" "$W/reads/run1/sa_1.fq" > fix/se/c.map
printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tPREFIX\ny\t%s\ty\n' "$W/fix/se/d" "$W/reads/run1/sb_1.fq" > fix/se/d.map
$MAPU merge --map fix/se/c.map fix/se/d.map --out "$W/fix/se/m2" > fix/se/merged2.map; note "    merge without SECOND: exit $?"

note "== 11. sample IDs that cannot name a file (map)"
{
  printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\n' "$W/fix/bad"
  i=0
  for id in 'a/b' 'c:d' '-e' '..' "$(printf 'f\x01g')" "$(printf 'h\xffi')" "$(head -c 201 /dev/zero | tr '\0' x)" 'ok_one'; do
    i=$((i + 1)); printf '%s\tgone_1.fq\tgone_2.fq\t%s\tp%s\n' "$id" "$W/run1/alignments/sa.sam.zst" "$i"
  done
} > fix/bad.map
p bad_sample_ids --map fix/bad.map --db "$DB" -t 4 --no_qcmsa
show bad_sample_ids "Line|Rename|Failed" 12
p bad_prefix --db "$DB" -t 4 --no_qcmsa --profile_only run1/alignments/sa.sam.zst,run1/alignments/sb.sam.zst --prefix 'a:b,-c' -o fix/bad_prefix
show bad_prefix "Error"
note "done"
