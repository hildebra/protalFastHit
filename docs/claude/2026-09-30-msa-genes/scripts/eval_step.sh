#!/bin/bash
# Evaluate a protal build + qcmsa on the strain audit's data, reusing its SAMs.
# Usage: eval_step.sh LABEL [acc] [phylo]   (binaries: ~/audit5/bin/protal_LABEL, ~/audit5/bin/qcmsa_LABEL.py)
#   acc:   accuracy runs A, B, Cs, Cl -> ~/audit5/accuracy/prot_<run>_LABEL, results/<run>_LABEL.*
#   phylo: known-tree runs -> ~/audit5/phylo/runs/<run>/protal_LABEL[k04], results_LABEL[k04].jsonl
set -u
L=$1; shift
parts=${*:-acc phylo}
BIN=$HOME/audit5/bin/protal_$L
QC=$HOME/audit5/bin/qcmsa_$L.py
DB=$HOME/audit5/world/protal_db
EXTRA=${EXTRA:-}
export PROTAL_QCMSA_SCRIPT=$QC

reprofile() {  # meta out_dir sam_dir log [extra args]
  local meta=$1 out=$2 sams=$3 log=$4; shift 4
  rm -rf "$out"; mkdir -p "$out/alignments"
  for f in "$sams"/*.sam.gz; do ln -s "$f" "$out/alignments/"; done
  $BIN --db $DB --map "$meta" -o "$out" -t 2 $EXTRA "$@" > "$log" 2>&1 || echo "protal $out exit $?"
}

case " $parts " in *" acc "*)
  A=$HOME/audit5/accuracy
  # The column-mapping copy of qcmsa must be the qcmsa protal ran.
  sed "s|~/audit5/src/scripts/qcmsa.py|$QC|" $A/make_qcmsa_colmap.py > $A/make_qcmsa_colmap_$L.py
  python3 $A/make_qcmsa_colmap_$L.py
  for run in A B Cs Cl; do
    man=$run; [ $run = Cs ] || [ $run = Cl ] && man=C
    out=$A/prot_${run}_$L; log=$A/logs/protal_${run}_$L.log
    reprofile $A/sim_$run/protal.meta $out $A/prot_$run/alignments $log
    python3 $A/evaluate.py $man $out ${run}_$L $log > $A/logs/eval_${run}_$L.log 2>&1 || echo "evaluate $run failed"
    echo "acc $run done"
  done
esac

case " $parts " in *" phylo "*)
  P=$HOME/audit5/phylo
  S2=$P/scripts_$L
  rm -rf $S2; cp -r $P/scripts $S2
  sed -i "s|~/audit5/src/scripts/qcmsa.py|$QC|" $S2/make_qcmsa_trace.py
  python3 $S2/make_qcmsa_trace.py > /dev/null
  one() {
    r=$1; D=$P/runs/$r; mix=""; sp="Malpha,Tone"
    [ $r = mixed ] && mix="--mix mix=t03"
    { [ $r = congener ] || [ $r = congener5 ]; } && sp="Malpha"
    for k in "" k04; do
      [ -z "$k" ] || [ -d $D/protal_k04 ] || continue
      [ -z "$k" ] || [ $sp != Malpha ] || continue
      sub=protal_$L$k
      sed "s#^\#OUTPUT_DIR\t.*#\#OUTPUT_DIR\t$D/$sub#" $D/sim/protal.meta > $D/$sub.meta
      if [ -z "$k" ]; then reprofile $D/$sub.meta $D/$sub $D/protal/alignments $D/$sub.log
      else reprofile $D/$sub.meta $D/$sub $D/protal/alignments $D/$sub.log --knob 0.4; fi
      if [ -z "$k" ]; then species=$sp; else species=Cferv; fi
      rm -rf $D/an_$L$k $D/results_$L$k.jsonl
      python3 $S2/analyze.py $D --species $species $mix --protal_subdir $sub --tag $L$k --variants raw,filt >> $D/eval_$L.log 2>&1
    done
    echo "phylo $r done"
  }
  export -f one reprofile; export P S2 L BIN DB EXTRA
  printf '%s\n' ${RUNS:-base20 d2 d2r2 d2r3 d3 d3r2 d5 d10 d50 uneven mixed congener congener5} | xargs -P 3 -I{} bash -c 'one {}'
esac
echo "EVAL $L DONE"
