#!/bin/bash
# The variants per read type: all features; without the ancestry sites; without the polymorphic group (with it the
# fixed-site features); without the column weights; without all three (every feature of the ancestral states);
# without the strain alleles. usage: ablate_all.sh <build outdir>
B=$1
S=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+gaps+untried
A=$(dirname $0)/ablate.sh
for t in pe pb; do
  bash $A $B $t full $S+ancestry+alleles+polymorphic+weights
  bash $A $B $t no_ancestry $S+alleles+polymorphic+weights
  bash $A $B $t no_polymorphic $S+ancestry+alleles+weights
  bash $A $B $t no_weights $S+ancestry+alleles+polymorphic
  bash $A $B $t no_ancestral $S+alleles
  bash $A $B $t no_alleles $S+ancestry+polymorphic+weights
done
echo ALL DONE
