#!/bin/bash
# The four variants per read type: all features, without the column weights, without the ancestry sites, without both.
B=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood
A=$(dirname $0)/ablate.sh
for t in pe pb se ont; do
  bash $A $t full $B+ancestry+gaps+untried+alleles+polymorphic+weights
  bash $A $t no_weights $B+ancestry+gaps+untried+alleles+polymorphic
  bash $A $t no_ancestry $B+gaps+untried+alleles+polymorphic+weights
  bash $A $t no_both $B+gaps+untried+alleles+polymorphic
done
echo ALL DONE
