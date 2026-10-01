#!/bin/bash
# Instructions of the alignment stage and the profiler with and without the gene neighbours (README.md), counted by
# callgrind on one thread: 40,000 read pairs, 400 PacBio and 600 Nanopore reads of the experiment's samples.
#   bash callgrind.sh WORK BUILD MODELS        (WORK: run.sh's folder, after it)
set -euo pipefail
W=$1 B=$2 M=$3
cd $W
head -n 160000 pe_R1.fq > cg_R1.fq
head -n 160000 pe_R2.fq > cg_R2.fq
head -n 1600 pb.fq > cg_pb.fq
head -n 2400 ont.fq > cg_ont.fq
models=(--model $M/trained_model.xml --model_pb $M/trained_model_pb.xml --model_ont $M/trained_model_ont.xml)
for type in pe pb ont; do
    reads=(-1 cg_$type.fq --read_type $type)
    stage='*RunLongReads*'
    [ $type = pe ] && reads=(-1 cg_R1.fq -2 cg_R2.fq) && stage='*RunPairedEnd*'
    for arm in with without; do
        extra=()
        [ $arm = without ] && extra=(--no_gene_neighbours)
        rm -rf cg_out_${type}_$arm
        # Everything is counted: the alignment runs in the OpenMP region's outlined function (*._omp_fn.*), which a
        # --toggle-collect on the run function does not reach; its inclusive cost is the alignment stage's.
        nice -n 10 valgrind --tool=callgrind --callgrind-out-file=cg_${type}_$arm.out $B/protal --db db "${reads[@]}" \
            --prefix $type -o cg_out_${type}_$arm -t 1 --no_strains --force "${models[@]}" "${extra[@]}" > cg_${type}_$arm.log 2>&1
        callgrind_annotate --inclusive=yes cg_${type}_$arm.out 2> /dev/null > cg_${type}_$arm.txt
        echo "$type $arm total $(grep -m1 'PROGRAM TOTALS' cg_${type}_$arm.txt | awk '{print $1}')" \
             "alignment $(grep -E "${stage//\*/}.*_omp_fn" cg_${type}_$arm.txt | head -1 | awk '{print $1}')" \
             "profiling $(grep -E 'ProfileSam' cg_${type}_$arm.txt | head -1 | awk '{print $1}')"
    done
done
