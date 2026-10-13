#!/bin/bash
B=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood
A=$(dirname $0)/ablate.sh
until grep -q "ALL DONE" $1; do sleep 20; done
for t in pb ont se; do bash $A $t no_ancestry $B+gaps+untried+alleles+polymorphic; done
echo MORE DONE
