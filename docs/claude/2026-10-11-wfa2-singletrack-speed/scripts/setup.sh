#!/usr/bin/env bash
# setup.sh: ~/wfasimd/src = git archive of audit-fixes HEAD (src.tar, made by the caller) with the scratch dump header
# (WFA2Wrapper2_dump.h: every WFA2 call written to $PROTAL_WFA_DUMP); protal and simulate_metagenomes built with -g;
# samples from the realworld_s1 genomes (pe 1M pairs, PacBio and Nanopore 30 Mb, 40 genomes); one protal run per
# sample (-t 1, --no_profile) that dumps its WFA2 calls to ~/wfasimd/calls/<sample>.tsv. On cores 0-3.
set -uo pipefail
S=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/6a0fc5a0-51ae-4ff5-a1a7-da994d3168e8/scratchpad/wfasimd
W=$HOME/wfasimd; E=$W/src; D=$W/data; C=$W/calls
DB=$HOME/realworld_s1/protal_db; TABLE=$HOME/realworld_s1/work/genomes_simulated.tsv
ONT_MODEL=$HOME/micromamba/envs/protal-db-build/data/QSHMM-ONT-HQ.model
pin() { taskset -c 0-3 nice -n 5 "$@"; }
mkdir -p $D $C
rm -rf $E; mkdir -p $E; tar -xf $W/src.tar -C $E
cp $S/WFA2Wrapper2_dump.h $E/src/Alignment/WFA2Wrapper2.h
cd $E
pin cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS=-g -DCMAKE_CXX_FLAGS=-g > $W/configure.log 2>&1 || { echo FAIL configure; exit 1; }
pin cmake --build build -j 4 --target protal simulate_metagenomes > $W/build.log 2>&1 && echo "OK build" || { echo FAIL build; grep -E ' error' -A3 $W/build.log | head -30; exit 1; }
SIM=$E/build/simulate_metagenomes
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
for s in pe pb ont; do
  case $s in
    pe) args="-1 $D/pe/reads/pe_s_1_R1.fq.gz -2 $D/pe/reads/pe_s_1_R2.fq.gz --read_type pe" ;;
    pb) args="-1 $D/pb.fq.gz --read_type pb" ;;
    ont) args="-1 $D/ont.fq.gz --read_type ont" ;;
  esac
  rm -rf $W/run_$s
  PROTAL_WFA_DUMP=$C/$s.tsv pin $E/build/protal --db $DB $args --prefix s -o $W/run_$s -t 1 --no_qcmsa --no_profile > $W/run_$s.log 2>&1
  echo "$s: exit $?, $(wc -l < $C/$s.tsv) WFA2 calls, $(du -h $C/$s.tsv | cut -f1)"
done
echo "SETUPDONE $(date +%T)"
