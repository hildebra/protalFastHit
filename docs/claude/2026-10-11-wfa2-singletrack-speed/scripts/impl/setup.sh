#!/usr/bin/env bash
# setup.sh: ~/wfaspeed/base = 0a1beaa (base.tar) with the dump header (WFA2Wrapper2_dump.h: every WFA2 call written to
# $PROTAL_WFA_DUMP, nothing else changed); ~/wfaspeed/new = the wfa-speed worktree's files (new.tar). Both Release -g;
# new with the unit tests, which run. Samples from the realworld_s1 genomes (pe 1M pairs, PacBio and Nanopore 30 Mb);
# one base run per sample (-t 1) that dumps its WFA2 calls to ~/wfaspeed/calls/<sample>.tsv and keeps its outputs for
# the identity check. On cores 0-3.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfaspeed
W=$HOME/wfaspeed; B=$W/base; N=$W/new; D=$W/data; C=$W/calls; O=$W/out
DB=$HOME/realworld_s1/protal_db; TABLE=$HOME/realworld_s1/work/genomes_simulated.tsv
ONT_MODEL=$HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model
pin() { taskset -c 0-3 nice -n 5 "$@"; }
mkdir -p $W $D $C $O
echo "== build base ($(date +%T))"
rm -rf $B; mkdir -p $B; tar -xf $S/base.tar -C $B
cp $S/WFA2Wrapper2_dump.h $B/src/Alignment/WFA2Wrapper2.h
pin cmake -S $B -B $B/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS=-g -DCMAKE_CXX_FLAGS=-g > $W/base.configure.log 2>&1 || { echo FAIL configure base; exit 1; }
pin cmake --build $B/build -j 4 --target protal simulate_metagenomes > $W/base.build.log 2>&1 && echo "OK base build" || { echo FAIL base build; grep -E ' error' -A3 $W/base.build.log | head -30; exit 1; }
echo "== build new ($(date +%T))"
rm -rf $N; mkdir -p $N; tar -xf $S/new.tar -C $N
pin cmake -S $N -B $N/build -G Ninja -DCMAKE_BUILD_TYPE=Release -DPROTAL_BUILD_TESTS=ON -DCMAKE_C_FLAGS=-g -DCMAKE_CXX_FLAGS=-g > $W/new.configure.log 2>&1 || { echo FAIL configure new; exit 1; }
pin cmake --build $N/build -j 4 --target protal simulate_metagenomes protal_tests > $W/new.build.log 2>&1 && echo "OK new build" || { echo FAIL new build; grep -E ' error|warning: ' -A3 $W/new.build.log | head -40; exit 1; }
grep -E "warning:" $W/new.build.log | grep -E "wfa2-lib|WFA2Wrapper2|test_WFA2Wrapper" | sort -u | head -20
echo "== unit tests, new ($(date +%T))"
pin $N/build/tests/protal_tests --gtest_filter='WFA2Wrapper.*:WfaTraceback.*' > $W/unit.wfa.log 2>&1; echo "wfa tests exit $?"
grep -E "^\[       OK|FAILED|^\[  PASSED|SKIPPED" $W/unit.wfa.log | head -30
(cd $N && pin ctest --test-dir build -j4 > $W/unit.all.log 2>&1); echo "ctest exit $?"
tail -5 $W/unit.all.log
echo "== samples ($(date +%T))"
SIM=$B/build/simulate_metagenomes
pin $SIM --genome_table $TABLE -o $D/pe -n 1 --sample_prefix pe_s --total_read_pairs 1000000 --species_per_sample 40 \
  --read_length 150 --seed 7 -t 4 --protal_metafile $D/pe/protal.meta > $D/pe.log 2>&1; echo "pe sim exit $?"
python3 - $D/pe/manifest.tsv $D <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
for name, bases, seed in (("pb", 30000000, 101), ("ont", 30000000, 202)):
    with open(f"{sys.argv[2]}/samples_{name}.tsv", "w") as fh:
        fh.write(f"sample\tout\tbases\tseed\n{name}\t{sys.argv[2]}/{name}.fq.gz\t{bases}\t{seed}\n")
    with open(f"{sys.argv[2]}/genomes_{name}.tsv", "w") as fh:
        fh.write("sample\tgenome\tfasta\tweight\thost\n")
        for r in rows:
            fh.write(f"{name}\t{r['genome']}\t{r['fasta_path']}\t{float(r['relative_abundance']) * float(r['genome_length'])!r}\t0\n")
EOF
pin $SIM --long_samples $D/samples_pb.tsv --long_genomes $D/genomes_pb.tsv --long_setup hifi:15000:3000:3 -t 4 > $D/pb.log 2>&1; echo "pb sim exit $?"
pin $SIM --long_samples $D/samples_ont.tsv --long_genomes $D/genomes_ont.tsv --long_setup qshmm:QSHMM-ONT-HQ:8000:6000:0.97:39/24/36 \
  --long_model $ONT_MODEL -t 4 > $D/ont.log 2>&1; echo "ont sim exit $?"
echo "== base runs with the calls dumped ($(date +%T))"
for s in pe pb ont; do
  case $s in
    pe) args="-1 $D/pe/reads/pe_s_1_R1.fq.gz -2 $D/pe/reads/pe_s_1_R2.fq.gz --read_type pe" ;;
    pb) args="-1 $D/pb.fq.gz --read_type pb" ;;
    ont) args="-1 $D/ont.fq.gz --read_type ont" ;;
  esac
  rm -rf $O/$s.base
  PROTAL_WFA_DUMP=$C/$s.tsv pin /usr/bin/time -f "TIME %e s, max RSS %M KB" $B/build/protal --db $DB $args --prefix s -o $O/$s.base -t 1 --no_qcmsa > $O/$s.base.log 2>&1
  echo "$s: exit $?, $(grep -m1 '^TIME' $O/$s.base.log), $(wc -l < $C/$s.tsv) WFA2 calls, $(du -h $C/$s.tsv | cut -f1)"
done
echo "SETUPDONE $(date +%T)"
