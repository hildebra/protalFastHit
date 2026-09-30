#!/bin/bash
# The "Invalid after alignment" diagnostics of a protal log: Info, Record and Anchor lines in full.
f=$1
grep -E '^(Info:|Record:|Anchor:)' $f | cut -c1-${2:-3000}
grep -n -A3 -B1 'Invalid' $f | cut -c1-200 | head -20
