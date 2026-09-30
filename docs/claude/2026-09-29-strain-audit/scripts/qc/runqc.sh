# usage: runqc.sh <strains dir> <out dir> [extra qcmsa args...]
S=$1; O=$2; shift 2
mkdir -p $O
for msa in $S/*.raw.msa.fna; do
  sp=$(basename $msa .raw.msa.fna)
  echo "=== $sp"
  python3 ~/audit5/src/scripts/qcmsa.py $msa $S/$sp.raw.partition.txt $S/$sp.meta.tsv --prefix $O/$sp --reapply-hcov 1000 "$@" 2>&1 | grep -v "^  iter\|^Saved\|^Done\|^Params"
  echo "rc=${PIPESTATUS[0]}"
  [ -f $O/$sp.msa.fna ] && echo "out: $(grep -c '>' $O/$sp.msa.fna) seqs, $(awk 'NR==2' $O/$sp.msa.fna | wc -c) first-line chars"
done
