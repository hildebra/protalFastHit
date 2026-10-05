#!/usr/bin/env bash
# The unit tests' summary for build <name>, and the ANISum of the line that varied (taxon 176, gene 90) per pe run.
W=$HOME/det-order; NAME=${1:-work}
$W/$NAME/tree/build/tests/protal_tests 2>&1 | grep -E 'tests ran|PASSED|FAILED'
for d in base/pe1 base/pe2 base/pe3 base/pe4 work/pe1 work/pe2; do
  echo "$d: $(awk -F'\t' '$2 == 176 && $6 == "Gene90" { print $10 }' $W/runs/$d/s.profile.gene.log)"
done
