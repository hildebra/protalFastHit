#!/usr/bin/env bash
# The whole-block prefetch (KmerLookupSM::PrefetchValues, PackedBlock::span) on HEAD: git archive HEAD with the two
# changed files (p2/src/Hash), Release; unit tests, the mini database and the e2e tests; the same outputs as HEAD's build
# (~/seed/head, built with -g: the same code) on the archaea world's 1M pairs and the small real world's pe, se, PacBio
# and ONT samples; seedbench rebuilt (its P0 keeps the old prefetch). 4 pinned cores.
set -uo pipefail
A=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/8e87f9e0-3682-453c-a7a2-636d4ec2ef61/scratchpad
W=$HOME/p2; REF=$HOME/seed/head/build/protal; N=$W/build/protal
pin=(taskset -c 0-3 nice -n 5)
rm -rf "$W"; mkdir -p "$W/src"; tar -xf "$A/archives/p2_head.tar" -C "$W/src"
cp "$A/p2/src/Hash/KmerLookup.h" "$A/p2/src/Hash/Seedmap.h" "$W/src/src/Hash/"
start=$(date +%s)
"${pin[@]}" cmake -S "$W/src" -B "$W/build" -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON > "$W/configure.log" 2>&1 &&
"${pin[@]}" cmake --build "$W/build" -j 4 --target protal protal_tests simulate_metagenomes > "$W/build.log" 2>&1 || { echo "BUILD FAIL"; grep -E "error" -A3 "$W/build.log" | head -30; exit 1; }
echo "build OK in $(( $(date +%s) - start )) s"
grep -E "warning" "$W/build.log" | grep -E "KmerLookup.h|Seedmap.h" | sort -u | head
"${pin[@]}" "$(find "$W/build" -name protal_tests -type f | head -1)" --gtest_brief=1 2>&1 | tail -3

compare() {  # a b: every file but misc/ and logs; SAM records sorted
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
R=$W/runs; mkdir -p "$R"
echo "[$(date +%T)] the archaea world, 1M pairs, aligning only"
S=$HOME/seed/samples
for v in ref new; do
  bin=$REF; [ $v = new ] && bin=$N
  "${pin[@]}" "$bin" --db "$HOME/seed/db/database.protal" -1 "$S/pe1m/reads/pe1m_1_R1.fq.gz" -2 "$S/pe1m/reads/pe1m_1_R2.fq.gz" -o "$R/arc_$v" \
      -t 4 --no_profile > "$R/arc_$v.log" 2>&1 || echo "FAILED arc $v"
done
compare "$R/arc_ref" "$R/arc_new"
echo "[$(date +%T)] the small real world: pe, se, PacBio, ONT, profiled"
SIM=$W/build/simulate_metagenomes; T=$W/samples; mkdir -p "$T/long"
"${pin[@]}" "$SIM" --genome_table "$HOME/realworld_s1/work/genomes.tsv" -o "$T/pe" -n 1 --sample_prefix pe --total_read_pairs 200000 \
    --species_per_sample 40 --strains_per_species 0.4 --congener_groups 0.5:2-5 --read_length 150 --sequencer HSXt --seed 5 -t 4 > "$T/pe.log" 2>&1
M=$T/pe/manifests/pe_1.tsv
for kind in pb ont; do
  printf 'sample\tout\tbases\tseed\n%s\t%s\t%s\t%s\n' "$kind" "$T/long/$kind.fq.gz" 20000000 9 > "$T/long/${kind}_samples.tsv"
  { printf 'sample\tgenome\tfasta\tweight\thost\n'; awk -F'\t' -v s="$kind" 'NR > 1 { print s "\t" $2 "\t" $11 "\t" $6 "\t0" }' "$M"; } > "$T/long/${kind}_genomes.tsv"
done
"${pin[@]}" "$SIM" --long_samples "$T/long/pb_samples.tsv" --long_genomes "$T/long/pb_genomes.tsv" --long_setup hifi:15000:3000:3 -t 4 > "$T/pb.log" 2>&1
"${pin[@]}" "$SIM" --long_samples "$T/long/ont_samples.tsv" --long_genomes "$T/long/ont_genomes.tsv" \
    --long_setup qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36 --long_model "$HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model" -t 4 > "$T/ont.log" 2>&1
DB=$HOME/realworld_s1/protal_db/database.protal
for s in pe se pb ont; do
  case $s in pe) reads=(-1 "$T/pe/reads/pe_1_R1.fq.gz" -2 "$T/pe/reads/pe_1_R2.fq.gz") ;; se) reads=(-1 "$T/pe/reads/pe_1_R1.fq.gz" --read_type se) ;;
             pb) reads=(-1 "$T/long/pb.fq.gz" --read_type pb) ;; ont) reads=(-1 "$T/long/ont.fq.gz" --read_type ont) ;; esac
  for v in ref new; do
    bin=$REF; [ $v = new ] && bin=$N
    "${pin[@]}" "$bin" --db "$DB" "${reads[@]}" -o "$R/${s}_$v" -t 4 --no_qcmsa > "$R/${s}_$v.log" 2>&1 || echo "FAILED $s $v"
  done
  compare "$R/${s}_ref" "$R/${s}_new"
done
grep -h "seeding:" "$R/arc_new.log" | cut -c1-200
echo "[$(date +%T)] mini database and e2e"
export PATH=$HOME/micromamba/envs/protal-db-build/bin:$PATH
cd "$W/src" || exit 1
PROTAL=$N "${pin[@]}" bash scripts/mini_db/build_mini_db.sh "$W/mini_db" > "$W/mini_db.log" 2>&1 || { echo "MINI DB FAIL"; tail "$W/mini_db.log"; }
PROTAL_TEST_DB=$W/mini_db/protal_db PROTAL=$N SIMULATE=$SIM "${pin[@]}" python3 -m unittest tests/e2e/test_protal_e2e.py > "$W/e2e.log" 2>&1
echo "e2e exit $?"; tail -3 "$W/e2e.log"; grep -a -E "^(FAIL|ERROR):" "$W/e2e.log" | head
echo "[$(date +%T)] seedbench against the new headers"
"${pin[@]}" bash "$A/seedbench_build.sh" "$A/seedbench.cpp" "$W/src" "$W/build" 4 > "$W/seedbench_build.log" 2>&1 && echo "seedbench built" || { echo "seedbench build FAIL"; tail "$W/seedbench_build.log"; }
echo "P2 CHECK DONE"
