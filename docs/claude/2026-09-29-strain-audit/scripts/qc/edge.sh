# Edge cases on copies of the lowcov run (its SAMs are reused, nothing is re-aligned).
set -u
Q=~/audit5/qc
SRC=$Q/runs/lowcov
P=~/audit5/bin/protal
export PROTAL_QCMSA_SCRIPT=$HOME/audit5/src/scripts/qcmsa.py
mk() {  # mk <name> : copy of the lowcov output dir + map pointing at it
  rm -rf $Q/runs/$1; mkdir -p $Q/runs/$1
  cp -r $SRC/protal $Q/runs/$1/protal
  sed "s#^\#OUTPUT_DIR.*#\#OUTPUT_DIR\t$Q/runs/$1/protal#" $SRC/sim/protal.meta > $Q/runs/$1/map.tsv
}
runp() { # runp <name> <label> args...
  n=$1; l=$2; shift 2
  $P --db ~/audit5/world/protal_db --map $Q/runs/$n/map.tsv -t 2 "$@" > $Q/runs/$n/$l.log 2>&1
  echo "[$n/$l] protal rc=$?"
  tr '\r' '\n' < $Q/runs/$n/$l.log | grep -a "qcmsa\]\|WARNING\|ERROR\|Error\|No gene\|No good\|across samples\|qcmsa.py:" | grep -v "^\[qcmsa\] '/" | head -20
}

echo "== E1 rerun with identical settings"
mk e1
md5sum $Q/runs/e1/protal/strains/*.fna | sort > $Q/runs/e1/before.md5
runp e1 rerun
md5sum $Q/runs/e1/protal/strains/*.fna | sort > $Q/runs/e1/after.md5
diff $Q/runs/e1/before.md5 $Q/runs/e1/after.md5 && echo "rerun: all MSAs identical"

echo "== E2 --msa_species naming a species that passes in no sample (--map_range 1-1: sample_1 lacks Mockella alpha)"
mk e2
touch -d '2020-01-01' $Q/runs/e2/protal/strains/s__Mockella_alpha.*
runp e2 nowhere --map_range 1-1 --msa_species s__Mockella_alpha
ls -l --time-style=+%Y-%m-%d $Q/runs/e2/protal/strains/s__Mockella_alpha.* | awk '{print $6, $5, $7}'
echo "meta lines: $(wc -l < $Q/runs/e2/protal/strains/s__Mockella_alpha.meta.tsv)"

echo "== E3 --msa_species naming a species that passes in one sample (--map_range 2-2)"
mk e3
runp e3 one --map_range 2-2 --msa_species s__Mockella_alpha
ls $Q/runs/e3/protal/strains/ | grep Mockella_alpha
grep '>' $Q/runs/e3/protal/strains/s__Mockella_alpha.raw.msa.fna

echo "== E4 qcmsa fails (bad flag)"
mk e4
runp e4 badflag --msa_species s__Mockella_alpha --qcmsa_args "--no-such-flag"
ls $Q/runs/e4/protal/strains/ | grep Mockella_alpha

echo "== E5 two map rows with the same #SAMPLEID"
mk e5
awk -F'\t' 'BEGIN{OFS="\t"} $1=="sample_3"{$1="sample_2"} {print}' $Q/runs/e5/map.tsv > $Q/runs/e5/m && mv $Q/runs/e5/m $Q/runs/e5/map.tsv
runp e5 dupid --msa_species s__Mockella_alpha
grep '>' $Q/runs/e5/protal/strains/s__Mockella_alpha.raw.msa.fna | tr '\n' ' '; echo
[ -f $Q/runs/e5/protal/strains/s__Mockella_alpha.msa.fna ] && grep '>' $Q/runs/e5/protal/strains/s__Mockella_alpha.msa.fna | tr '\n' ' '; echo
iqtree2 -s $Q/runs/e5/protal/strains/s__Mockella_alpha.msa.fna -m GTR+G -B 1000 -T 2 --seqtype DNA --prefix $Q/runs/e5/iq -redo 2>&1 | grep -i "error\|duplicat" | head -3

echo "== E6 #SAMPLEID with a space and parentheses"
mk e6
awk -F'\t' 'BEGIN{OFS="\t"} $1=="sample_3"{$1="pt 3 (day 1)"} {print}' $Q/runs/e6/map.tsv > $Q/runs/e6/m && mv $Q/runs/e6/m $Q/runs/e6/map.tsv
runp e6 oddid --msa_species s__Mockella_alpha
grep '>' $Q/runs/e6/protal/strains/s__Mockella_alpha.msa.fna | tr '\n' '|'; echo
iqtree2 -s $Q/runs/e6/protal/strains/s__Mockella_alpha.msa.fna -m GTR+G -B 1000 -T 2 --seqtype DNA --prefix $Q/runs/e6/iq -redo 2>&1 | grep -i -A3 "error\|changed\|illegal" | head -8
cat $Q/runs/e6/iq.treefile 2>/dev/null
