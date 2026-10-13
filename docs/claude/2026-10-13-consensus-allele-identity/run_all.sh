#!/bin/bash
# The comparison: the reduced database built four times from the same inputs and seed, and the mock samples profiled
# with each.
#   base       protal as committed (49cd583)
#   both       the consensus tie-break and the allele identities (models trained on the default set, with them)
#   consensus  the same protal, models trained without the allele identities (--features: the set before them)
#   identity   the allele identities alone (protal without the tie-break, the default set)
# Needs the clones ~/idtest/{base,new,idonly} (setup.sh) and the inputs ~/idtest/inputs
# (../2026-10-13-mock-local-test/make_inputs.sh). usage: run_all.sh     (WSL; 4 cores)
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
W=$HOME/idtest
L=$W/logs
BEFORE=normalized+adjacency+distance+depth+divergence+unfiltered+ref+complexity+consistency+shape+neighbourhood+ancestry+gaps+untried+alleles+polymorphic+weights
for spec in both:new base:base identity:idonly consensus:new; do
  IFS=: read -r name clone <<< "$spec"
  extra=""; [ $name = consensus ] && extra="--features $BEFORE"
  if [ ! -f $W/db_$name/protal_db/database.protal ]; then
    echo "build $name $(date +%T)"
    bash $HERE/build_db.sh $W/$clone $W/inputs $W/db_$name $extra > $L/build_$name.log 2>&1 || echo "build $name failed: $?"
  fi
  echo "mocks $name $(date +%T)"
  bash $HERE/run_mocks.sh $W/$clone $W/db_$name/protal_db idc_$name > $L/mocks_$name.log 2>&1 || echo "mocks $name failed: $?"
done
echo "all done $(date +%T)"
