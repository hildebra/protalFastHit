#!/usr/bin/env bash
# Tables read only where the run reads them: HEAD (0a1beaa, ~/vc from verify_commit.sh) against HEAD + the change
# (~/tl/work), on the small real world's database (~/realworld_s1: strain alleles, column weights, suspect copies). Same
# outputs for pe (alleles before the alignment, as before), pe --no_allele_scores and PacBio (alleles and weights beside
# it now), a PacBio + ONT cohort with --profile_ahead, PacBio --no_profile (the SAM), and a profile-only run with
# --keep_suspect_copies; then the unit tests and the e2e tests on a fresh mini database. 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
REPO=/mnt/c/Users/hildebra/Documents/locDev/protal
until grep -q "VERIFY DONE" "$A/verify_commit.out" 2> /dev/null; do sleep 15; done
T=$HOME/tl; REF=$HOME/vc/build/protal; N=$T/work/build/protal
pin=(taskset -c 0-3 nice -n 5)
rm -rf "$T"; mkdir -p "$T/work/src"; tar -xf "$A/archives/commit.tar" -C "$T/work/src"
for f in src/RunProtal.h tests/e2e/test_protal_e2e.py; do cp "$REPO/$f" "$T/work/src/$f"; done
start=$(date +%s)
"${pin[@]}" cmake -S "$T/work/src" -B "$T/work/build" -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > "$T/configure.log" 2>&1 &&
"${pin[@]}" cmake --build "$T/work/build" -j 4 --target protal protal_tests simulate_metagenomes > "$T/build.log" 2>&1 || { echo "BUILD FAIL"; grep -E "error" -A3 "$T/build.log" | head -30; exit 1; }
echo "build OK in $(( $(date +%s) - start )) s"
grep -E "warning" "$T/build.log" | grep -E "RunProtal.h" | sort -u | head
"${pin[@]}" "$(find "$T/work/build" -name protal_tests -type f | head -1)" --gtest_brief=1 2>&1 | tail -3

# Samples from the small real world's genomes.
S=$T/samples; SIM=$T/work/build/simulate_metagenomes; mkdir -p "$S/long"
"${pin[@]}" "$SIM" --genome_table "$HOME/realworld_s1/work/genomes.tsv" -o "$S/pe" -n 1 --sample_prefix pe --total_read_pairs 200000 \
    --species_per_sample 40 --strains_per_species 0.4 --congener_groups 0.5:2-5 --read_length 150 --sequencer HSXt --seed 5 -t 4 > "$S/pe.log" 2>&1
M=$S/pe/manifests/pe_1.tsv
for kind in pb ont; do
  printf 'sample\tout\tbases\tseed\n%s\t%s\t%s\t%s\n' "$kind" "$S/long/$kind.fq.gz" 20000000 9 > "$S/long/${kind}_samples.tsv"
  { printf 'sample\tgenome\tfasta\tweight\thost\n'; awk -F'\t' -v s="$kind" 'NR > 1 { print s "\t" $2 "\t" $11 "\t" $6 "\t0" }' "$M"; } > "$S/long/${kind}_genomes.tsv"
done
"${pin[@]}" "$SIM" --long_samples "$S/long/pb_samples.tsv" --long_genomes "$S/long/pb_genomes.tsv" --long_setup hifi:15000:3000:3 -t 4 > "$S/pb.log" 2>&1
"${pin[@]}" "$SIM" --long_samples "$S/long/ont_samples.tsv" --long_genomes "$S/long/ont_genomes.tsv" \
    --long_setup qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36 --long_model "$HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model" -t 4 > "$S/ont.log" 2>&1
pe1=$S/pe/reads/pe_1_R1.fq.gz pe2=$S/pe/reads/pe_1_R2.fq.gz
DB=$HOME/realworld_s1/protal_db/database.protal
R=$T/runs; mkdir -p "$R"
run() {  # name binary args...
  local name=$1 bin=$2; shift 2
  rm -rf "$R/$name"
  "${pin[@]}" "$bin" --db "$DB" -o "$R/$name" -t 4 --no_qcmsa "$@" > "$R/$name.log" 2>&1 || echo "FAILED $name"
}
cohort() {  # name binary opts...
  local name=$1 bin=$2; shift 2
  rm -rf "$R/$name"; mkdir -p "$R/$name"
  printf '#OUTPUT_DIR\t%s\n#SAMPLEID\tFIRST\tSECOND\tPREFIX\tREAD_TYPE\npb\t%s\t-\tpb\tpb\nont\t%s\t-\tont\tont\n' "$R/$name" "$S/long/pb.fq.gz" "$S/long/ont.fq.gz" > "$R/$name.map"
  "${pin[@]}" "$bin" --db "$DB" --map "$R/$name.map" -t 4 --no_qcmsa "$@" > "$R/$name.log" 2>&1 || echo "FAILED $name"
}
for v in ref work; do
  bin=$REF; [ $v = work ] && bin=$N
  run "${v}_pe" "$bin" -1 "$pe1" -2 "$pe2"
  run "${v}_pe_noscores" "$bin" -1 "$pe1" -2 "$pe2" --no_allele_scores
  run "${v}_pb" "$bin" -1 "$S/long/pb.fq.gz" --read_type pb
  run "${v}_pb_noprofile" "$bin" -1 "$S/long/pb.fq.gz" --read_type pb --no_profile
  cohort "${v}_long_ahead" "$bin" --profile_ahead
done
sam=$(find "$R/ref_pe" -maxdepth 1 -name "*.sam.zst" | head -1)
for v in ref work; do
  bin=$REF; [ $v = work ] && bin=$N
  rm -rf "$R/${v}_keep"; "${pin[@]}" "$bin" --db "$DB" --profile_only "$sam" -o "$R/${v}_keep" -t 4 --no_qcmsa --keep_suspect_copies > "$R/${v}_keep.log" 2>&1 || echo "FAILED ${v}_keep"
done
compare() {  # a b
  local a=$1 b=$2 bad=0
  while IFS= read -r f; do
    rel=${f#$a/}
    case $rel in misc/*|*.log) continue ;; esac
    if [[ $rel == *.sam.zst ]]; then
      cmp -s <(zstd -dc "$a/$rel" | grep -v '^@' | sort) <(zstd -dc "$b/$rel" 2> /dev/null | grep -v '^@' | sort) || { echo "  differs: $rel"; bad=1; }
    else
      cmp -s "$a/$rel" "$b/$rel" || { echo "  differs: $rel"; bad=1; }
    fi
  done < <(find "$a" -type f | sort)
  echo "$([ $bad = 0 ] && echo same || echo DIFFERENT): $(basename "$a") vs $(basename "$b") ($(find "$a" -type f | grep -v '/misc/' | wc -l) files)"
}
for x in pe pe_noscores pb pb_noprofile long_ahead keep; do compare "$R/ref_$x" "$R/work_$x"; done
for x in pe pe_noscores pb pb_noprofile long_ahead keep; do
  echo "== work_$x"
  grep -E "^(Strain alleles|Column weights|Suspect copies):|^Loading the (taxonomy|profiling)" "$R/work_$x.log" | cut -c1-170
done
echo "[$(date +%T)] e2e"
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
cd "$T/work/src" || exit 1
PROTAL=$N "${pin[@]}" bash scripts/mini_db/build_mini_db.sh "$T/mini_db" > "$T/mini_db.log" 2>&1 || { echo "MINI DB FAIL"; tail "$T/mini_db.log"; }
PROTAL_TEST_DB=$T/mini_db/protal_db PROTAL=$N SIMULATE=$T/work/build/simulate_metagenomes "${pin[@]}" python3 -m unittest tests/e2e/test_protal_e2e.py > "$T/e2e.log" 2>&1
echo "e2e exit $?"; tail -3 "$T/e2e.log"; grep -a -E "^(FAIL|ERROR):" "$T/e2e.log" | head
echo "TL DONE"
