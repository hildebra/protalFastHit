#!/bin/bash
# Two more variants for pe and pb, for the cost against the gain (section 5): without every feature that needs the
# ancestry sites (ancestry, polymorphic's ancestry_fixed_*, weights' ancestry_agreement_weighted: the three groups
# whole, as the trainer takes groups), and without the two groups that cost most (ancestry, distance).
A=${ABLATE:-$(dirname $0)/ablate.sh}
FGROUPS="normalized adjacency distance depth divergence unfiltered ref complexity consistency shape neighbourhood ancestry gaps untried alleles polymorphic weights"
without() {
  local out=""
  for g in $FGROUPS; do
    skip=0; for x in "$@"; do [ $g = $x ] && skip=1; done
    [ $skip = 0 ] && out="$out+$g"
  done
  echo ${out#+}
}
for t in ${TYPES:-pe pb}; do
  bash $A $t no_sites $(without ancestry polymorphic weights)
  bash $A $t no_anc_dist $(without ancestry distance)
done
echo EXTRA DONE
