#!/bin/bash
# The experiment of README.md: an operon-like world, a database with gene neighbours, paired-end, PacBio and
# Nanopore samples of it, and protal on each with and without --no_gene_neighbours (3 alternated runs each).
#   bash run.sh WORK SRC BUILD MODELS
# WORK: the experiment's folder; SRC: a protal checkout; BUILD: its build folder (protal); MODELS: a folder with
# trained_model.xml, trained_model_pb.xml, trained_model_ont.xml (the V3 tuning build's). Needs pbsim (pbsim3) and
# its models (data/ or share/pbsim*/ next to its bin/). Runs niced: the machine is shared.
set -euo pipefail
W=$1 SRC=$2 B=$3 M=$4
S=$SRC/scripts/mini_db
here=$(cd "$(dirname "$0")" && pwd)
T=4
mkdir -p $W
cd $W

if [ ! -s gtdb/simulation/genomes.tsv ]; then
    python3 $S/gtdb_like_lineages.py --species 80 --archaea 0.1 --seed 3 > lineages.txt
    python3 $S/simulate_gtdb_release.py --outdir gtdb --lineages lineages.txt --operons --genome_length 400000 \
        --genomes_per_species 2 --contigs 5 --strain_divergence 0.005-0.03 --species_divergence 0.015-0.06 \
        --gene_rates categories --seed 3 2> simulate.log
fi
if [ ! -s db/database.protal ]; then
    rm -rf conv db
    python3 $S/gtdb_to_protal_db.py --gtdb gtdb --outdir conv -t $T > convert.log 2>&1
    python3 $S/gene_neighbours.py --db conv --genome_table gtdb/simulation/genomes.tsv -t $T --positions positions.tsv > gene_neighbours.log
    python3 $here/design.py heldout gtdb/simulation/genomes.tsv > heldout.txt
    python3 $S/gtdb_to_protal_db.py --from_db conv --exclude_species heldout.txt --outdir db -t $T > exclude.log 2>&1
    full=$(ls db/full_reference.fna* | head -1)
    nice -n 10 $B/protal --build --no_profile -t $T --db db --reference db/reference.fna --full_reference $full > build.log 2>&1
fi

# The samples: 22 species of the database (half from their representative's genome, half from the other genome,
# a strain 0.5-3% from it) and the 8 held out, lognormal abundances.
if [ ! -s pe_R1.fq ]; then
    python3 $here/design.py community gtdb/simulation/genomes.tsv heldout.txt > community.tsv
    python3 $S/simulate_reads.py --genomes gtdb/simulation/genomes.tsv --community community.tsv --out_prefix pe \
        --pairs 300000 --seed 9 > reads.log 2>&1
    python3 $here/design.py truth community.tsv gtdb/simulation/genomes.tsv conv/internal_taxonomy.dmp > truth.tsv
fi
pbsim_data=$(dirname "$(dirname "$(command -v pbsim)")")  # the models are in its data/ or share/pbsim*/
long_reads() {  # type method model length sd accuracy
    local type=$1 method=$2 model=$3 length=$4 sd=$5 accuracy=$6
    [ -s $type.fq ] && return
    local modelfile
    modelfile=$(find $pbsim_data/data $pbsim_data/share -name "$model.model" 2>/dev/null | head -1)
    [ -n "$modelfile" ] || { echo "no pbsim3 model $model under $pbsim_data" >&2; exit 1; }
    rm -rf lr_$type && mkdir lr_$type
    local mean  # depth: 2 on average over the genomes, by abundance
    mean=$(tail -n +2 community.tsv | awk '{s += $2} END {print s / NR}')
    tail -n +2 community.tsv | while IFS=$'\t' read -r acc abundance; do
        fasta=$(awk -F'\t' -v a=$acc '$1 == a {print $3}' gtdb/simulation/genomes.tsv)
        zcat -f $fasta > lr_$type/genome.fna
        depth=$(python3 -c "print(round(max(0.3, 2 * $abundance / $mean), 3))")
        pbsim --strategy wgs --method $method --$method $modelfile --genome lr_$type/genome.fna --depth $depth \
            --length-mean $length --length-sd $sd --accuracy-mean $accuracy --seed 5 --prefix lr_$type/$acc \
            --id-prefix $acc >> lr_$type/pbsim.log 2>&1
    done
    zcat -f lr_$type/*.fq* lr_$type/*.fastq* 2>/dev/null > $type.fq || true  # pbsim3 writes .fq.gz or .fastq
    [ -s $type.fq ] || { echo "pbsim wrote no reads (lr_$type/pbsim.log)" >&2; exit 1; }
    rm -rf lr_$type
}
long_reads pb errhmm ERRHMM-SEQUEL 15000 3000 0.999
long_reads ont qshmm QSHMM-ONT-HQ 8000 6000 0.97

run() {  # name type with|without
    local name=$1 type=$2 arm=$3 extra=()
    [ $arm = without ] && extra=(--no_gene_neighbours)
    local reads=(-1 $type.fq --read_type $type)
    [ $type = pe ] && reads=(-1 pe_R1.fq -2 pe_R2.fq)
    rm -rf out/$name
    nice -n 10 /usr/bin/time -f "%e %U %M" -o time_$name.txt $B/protal --db db "${reads[@]}" --prefix $type -o out/$name \
        -t $T --no_strains --profile_truth truth.tsv --model $M/trained_model.xml --model_pb $M/trained_model_pb.xml \
        --model_ont $M/trained_model_ont.xml --force "${extra[@]}" > log_$name.txt 2>&1
}
for rep in 1 2 3; do
    for type in pe pb ont; do
        for arm in with without; do
            run ${type}_${arm}_$rep $type $arm
        done
    done
done
python3 $here/design.py report . > results.md
cat results.md
