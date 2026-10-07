#!/usr/bin/env bash
# run.sh - combining the .sam.zst files of several protal runs into one strain analysis
# (docs/claude/2026-10-07-sam-combine). Three studies are aligned and profiled on their own (maps from
# protal_map_utils generate), sa..sd once more in one joint run as the reference, and then every way there is to
# combine the studies' SAMs is tried. Every protal run on 4 cores (taskset -c 0-3), --no_qcmsa.
#
#   SRC=<git archive of the commit> PROTAL=<its protal> DB=<mini database.protal> W=<work dir> bash run.sh
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
SRC=${SRC:-$HOME/samcombine/head}
PROTAL=${PROTAL:-$HOME/samcombine/bin/protal}
DB=${DB:-$HOME/samcombine/db/database.protal}
W=${W:-$HOME/samcombine/work}
RUN="taskset -c 0-3 nice -n 5"
MAPU="python3 $SRC/scripts/protal_map_utils"

rm -rf "$W"; mkdir -p "$W/logs"; cd "$W" || exit 1
SUMMARY=$W/summary.txt; : > "$SUMMARY"

note() { echo "$*" | tee -a "$SUMMARY"; }
# p NAME ARGS...: protal in ${CWD:-$W}; its output in logs/NAME.log, its exit code in the summary
p() {
  local name=$1; shift
  (cd "${CWD:-$W}" && $RUN "$PROTAL" "$@") > "$W/logs/$name.log" 2>&1
  local rc=$?
  note "$name: exit $rc"
}
# last NAME [N]: the last lines of a run's log that say what went wrong, or its last N lines
last() { grep -E "Error|error|rror:|Warning|terminate|what\(\)|Abort|All alignments|Skip " "$W/logs/$1.log" | sort | uniq -c | head -${2:-12} | sed 's/^/    /' | tee -a "$SUMMARY"; }
# msa_rows DIR: each raw MSA's species and row names
msa_rows() {
  for f in "$1"/*.raw.msa.fna; do
    [ -e "$f" ] || { note "    (no raw MSAs in $1)"; return; }
    note "    $(basename "$f" .raw.msa.fna): $(grep '>' "$f" | tr -d '>' | tr '\n' ' ')"
  done
}
# same_msas A B: the raw MSAs, partitions and meta tables of strain folder A compared with B's
same_msas() {
  local diffs=0 n=0
  for f in "$1"/*.raw.msa.fna "$1"/*.raw.partition.txt "$1"/*.meta.tsv; do
    [ -e "$f" ] || continue
    n=$((n + 1))
    cmp -s "$f" "$2/$(basename "$f")" || { diffs=$((diffs + 1)); note "    differs: $(basename "$f")"; }
  done
  note "    $n strain files of $1 against $2: $diffs differ"
}

PROTAL_TEST_DB=$DB PROTAL=$PROTAL python3 "$HERE/simulate.py" reads --root "$SRC" > logs/simulate.log 2>&1 \
  || { cat logs/simulate.log; exit 1; }

note "== three studies, each run on its own; and sa..sd in one joint run"
for r in run1 run2 run3; do
  $MAPU generate --input "reads/$r" --out "$r" > "$r.map"
  p "$r" --map "$r.map" --db "$DB" -t 4 --no_qcmsa
done
$MAPU generate --input reads/run1 reads/run2 --out joint > joint.map
p joint --map joint.map --db "$DB" -t 4 --no_qcmsa
msa_rows joint/strains
note "    run1/alignments: $(ls run1/alignments | tr '\n' ' ')"

note "== A. --profile_only of run1's and run2's SAMs (a shell glob joined with commas), -o comb_po"
touch -d '2000-01-01' run1/alignments/sa.sam.zst.err
SAMS=$(ls run1/alignments/*.sam.zst run2/alignments/*.sam.zst | paste -sd,)
p A_profile_only --db "$DB" --profile_only "$SAMS" -o comb_po -t 4 --no_qcmsa
msa_rows comb_po/strains
same_msas joint/strains comb_po/strains
for s in sa sb sc sd; do
  cmp -s "joint/profiles/$s.profile" "comb_po/$s.profile" && continue
  note "    profile of $s differs from the joint run's:"
  diff "joint/profiles/$s.profile" "comb_po/$s.profile" | head -6 | sed 's/^/      /' | tee -a "$SUMMARY"
done
note "    run1/alignments/sa.sam.zst.err rewritten by A: $(stat -c %y run1/alignments/sa.sam.zst.err | cut -c1-19)"
note "    A2. the joint run's own SAMs, --profile_only:"
p A2_profile_only_joint --db "$DB" --profile_only "$(ls joint/alignments/*.sam.zst | paste -sd,)" -o comb_po_joint -t 4 --no_qcmsa
same_msas joint/strains comb_po_joint/strains
note "    A3. A again into the same -o: rows of comb_po/misc/sa_runtime.tsv before: $(wc -l < comb_po/misc/sa_runtime.tsv)"
p A3_profile_only_again --db "$DB" --profile_only "$SAMS" -o comb_po -t 4 --no_qcmsa
note "    after: $(wc -l < comb_po/misc/sa_runtime.tsv)"

note "== B. --profile_only without -o (run in an empty folder B)"
mkdir -p B
CWD=$W/B p B_no_outdir --db "$DB" --profile_only "$W/run1/alignments/sa.sam.zst,$W/run2/alignments/sc.sam.zst" -t 4 --no_qcmsa
last B_no_outdir
note "    B/: $(ls B | tr '\n' ' ')   B/strains: $(ls B/strains 2>/dev/null | wc -l) files"
note "    new in run1/alignments: $(ls run1/alignments | grep -v '\.sam\.zst' | tr '\n' ' ')"
tail -3 logs/B_no_outdir.log | sed 's/^/    | /' | tee -a "$SUMMARY"

note "== C. two studies' samples of the same name (run1/alignments/sa and run3/alignments/sa), -o comb_c"
p C_same_name --db "$DB" --profile_only "run1/alignments/sa.sam.zst,run3/alignments/sa.sam.zst,run2/alignments/sc.sam.zst" -o comb_c -t 4 --no_qcmsa
last C_same_name
p C_same_name_prefix --db "$DB" --profile_only "run1/alignments/sa.sam.zst,run3/alignments/sa.sam.zst,run2/alignments/sc.sam.zst" \
  --prefix run1_sa,run3_sa,sc -o comb_c2 -t 4 --no_qcmsa
msa_rows comb_c2/strains

note "== D. a folder, and a glob protal would have to expand itself"
p D_folder --db "$DB" --profile_only run1/alignments -o comb_d -t 4 --no_qcmsa
tail -2 logs/D_folder.log | sed 's/^/    | /' | tee -a "$SUMMARY"
p D_glob --db "$DB" --profile_only 'run*/alignments/*.sam.zst' -o comb_d2 -t 4 --no_qcmsa
last D_glob 4

note "== E. protal_map_utils merge of run1.map and run2.map"
$MAPU merge --map run1.map run2.map --out "$W/comb_merge" > merged.map; note "    merge: exit $?"
sed 's/^/    | /' merged.map | tee -a "$SUMMARY"
p E_merge --map merged.map --db "$DB" -t 4 --no_qcmsa
note "    'All alignments are present': $(grep -c 'All alignments are present' logs/E_merge.log); index loaded: $(grep -c 'Index loaded\|Load index\|Loading the index' logs/E_merge.log)"
note "    comb_merge/alignments: $(ls comb_merge/alignments 2>/dev/null | tr '\n' ' ')"
same_msas joint/strains comb_merge/strains
$MAPU merge --map run1.map run2.map > merged_noout.map; note "    merge without --out: exit $?"
sed 's/^/    | /' merged_noout.map | tee -a "$SUMMARY"
$MAPU merge --map run1.map run2.map --out "$W/comb_merge2" > merged2.map
mv reads reads.away
p E_merge_without_reads --map merged2.map --db "$DB" -t 4 --no_qcmsa
last E_merge_without_reads 6

note "== F. a map of the SAMs (absolute SAM paths, read files that are gone), #OUTPUT_DIR comb_map"
{
  printf '#OUTPUT_DIR\t%s\n' "$W/comb_map"
  printf '#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\n'
  for s in run1/alignments/sa run1/alignments/sb run2/alignments/sc run2/alignments/sd; do
    printf '%s\tgone_R1.fq\tgone_R2.fq\t%s\t%s\n' "$(basename "$s")" "$W/$s.sam.zst" "$(basename "$s")"
  done
} > sams.map
p F_map_of_sams --map sams.map --db "$DB" -t 4 --no_qcmsa
last F_map_of_sams 6
same_msas joint/strains comb_map/strains
note "    F2. with run3's sa as study3_sa (#SAMPLEID and PREFIX), the same file name as run1's sa"
{
  printf '#OUTPUT_DIR\t%s\n' "$W/comb_map2"
  printf '#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\n'
  printf 'study1_sa\t-\t-\t%s\tstudy1_sa\n' "$W/run1/alignments/sa.sam.zst"
  printf 'study3_sa\t-\t-\t%s\tstudy3_sa\n' "$W/run3/alignments/sa.sam.zst"
  printf 'study2_sc\t-\t-\t%s\tstudy2_sc\n' "$W/run2/alignments/sc.sam.zst"
} > sams2.map
p F2_map_same_names --map sams2.map --db "$DB" -t 4 --no_qcmsa
last F2_map_same_names 6
msa_rows comb_map2/strains
mv reads.away reads

note "== G. a map of two existing SAMs and one new sample's reads"
{
  printf '#OUTPUT_DIR\t%s\n' "$W/comb_g"
  printf '#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\n'
  printf 'sa\tgone_R1.fq\tgone_R2.fq\t%s\tsa\n' "$W/run1/alignments/sa.sam.zst"
  printf 'sb\tgone_R1.fq\tgone_R2.fq\t%s\tsb\n' "$W/run1/alignments/sb.sam.zst"
  printf 'sc\t%s\t%s\tsc.sam.zst\tsc\n' "$W/reads/run2/sc_1.fq" "$W/reads/run2/sc_2.fq"
} > grow.map
p G_add_a_sample --map grow.map --db "$DB" -t 4 --no_qcmsa
last G_add_a_sample 8
msa_rows comb_g/strains

note "== H. SAMs in a read-only folder (run2/alignments), --profile_only -o comb_h"
rm -f run2/alignments/*.err; chmod a-w run2/alignments
p H_read_only --db "$DB" --profile_only "$(ls run2/alignments/*.sam.zst | paste -sd,)" -o comb_h -t 4 --no_qcmsa
last H_read_only 4
note "    .err files in run2/alignments: $(ls run2/alignments | grep -c '\.err$'); comb_h: $(ls comb_h | tr '\n' ' ')"
chmod u+w run2/alignments

note "== L. --map with --profile_only (run1.map)"
p L_map_and_profile_only --map run1.map --profile_only ignored.sam.zst --db "$DB" -t 4 --no_qcmsa
last L_map_and_profile_only 6

note "== N. protal_map_utils generate: the sample IDs of common read file names (empty files)"
for names in "x_R1.fq x_R2.fq" "x_1.fq x_2.fq" "x.R1.fastq.gz x.R2.fastq.gz" "x_S1_L001_R1_001.fastq.gz x_S1_L001_R2_001.fastq.gz" \
             "x_1.fq.zst x_2.fq.zst"; do
  rm -rf names && mkdir names
  for f in $names; do : > "names/$f"; done
  note "    $names -> $($MAPU generate --input names 2>&1 | awk -F'\t' '!/^#/ {print $1}' | tr '\n' ' ')"
done

note "done"
