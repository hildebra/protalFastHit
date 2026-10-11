#!/bin/bash
# Every feature group of the default set left out in turn, per read type (pe first, then pb, se, ont), and three
# combinations: the four strain-level groups of 2026-10-08..10 (ancestry, alleles, polymorphic, weights) together,
# those plus gaps and untried, and the old core set alone (normalized..unfiltered).
A=$(dirname $0)/ablate.sh
FGROUPS="normalized adjacency distance depth divergence unfiltered ref complexity consistency shape neighbourhood ancestry gaps untried alleles polymorphic weights"
without() {  # the default set without the groups named
  local out=""
  for g in $FGROUPS; do
    skip=0; for x in "$@"; do [ $g = $x ] && skip=1; done
    [ $skip = 0 ] && out="$out+$g"
  done
  echo ${out#+}
}
for t in ${TYPES:-pe pb se ont}; do
  bash $A $t full $(without)
  for g in $FGROUPS; do [ $g = normalized ] || bash $A $t no_$g $(without $g); done
  bash $A $t no_strain4 $(without ancestry alleles polymorphic weights)
  bash $A $t no_strain6 $(without ancestry alleles polymorphic weights gaps untried)
  bash $A $t core normalized+adjacency+distance+depth+divergence+unfiltered
done
echo ALL DONE
