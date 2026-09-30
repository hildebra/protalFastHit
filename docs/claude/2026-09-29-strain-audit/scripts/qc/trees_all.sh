# usage: trees_all.sh <run name> <strains dir> <manifest> <ablation dir>
name=$1; S=$2; M=$3; A=$4
cd ~/audit5/qc
args=()
for msa in $S/*.raw.msa.fna; do
  sp=$(basename $msa .raw.msa.fna)
  args+=("raw=$msa")
  for c in default min-parsimony_0; do
    [ -f $A/$c/$sp.msa.fna ] && args+=("$c=$A/$c/$sp.msa.fna")
  done
done
python3 treeeval.py $M trees/$name "${args[@]}"
