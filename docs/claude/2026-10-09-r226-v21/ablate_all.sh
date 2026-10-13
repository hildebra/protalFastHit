#!/bin/bash
B=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood
A=$(dirname $0)/ablate.sh
for t in pe pb ont se; do
  bash $A $t full $B+ancestry+gaps+untried+alleles+polymorphic
  bash $A $t no_new $B+ancestry+gaps+untried
  bash $A $t no_ancestry_new $B+gaps+untried
  [ $t = pe ] && bash $A $t no_poly $B+ancestry+gaps+untried+alleles
  [ $t = pe ] && bash $A $t no_ancestry $B+gaps+untried+alleles+polymorphic
done
echo ALL DONE
