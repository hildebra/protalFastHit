#!/bin/bash
# The suites as CI runs them, old tree (HEAD) against new (the changes), alternated: old new, new old.
SP=/mnt/c/Users/hildebra/AppData/Local/Temp/claude/C--Users-hildebra-Documents-locDev-protal/65c3e396-1c71-480e-8173-2f2320bc22de/scratchpad
A=~/testspeed/ab; rm -rf $A; mkdir -p $A
old() { bash $SP/run_all.sh ~/testspeed/src ~/testspeed/build $A/old$1 unit mini e2e example minidb scripts; }
new() { bash $SP/run_all.sh ~/testspeed/new ~/testspeed/build $A/new$1 unit ciexample e2e minidb scripts; }
old 1; new 1; new 2; old 2
echo finished > $A/finished
