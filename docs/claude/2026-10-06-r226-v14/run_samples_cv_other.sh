#!/bin/bash
# samples_cv.py for se, pb and ont (v14 as built, and + the sample's complexity), each at its build's knob.
cd /mnt/c/Users/hildebra/Documents/locDev/protal/docs/claude/2026-10-06-r226-v14 || exit 1
for spec in se:0.73 pb:0.5 ont:0.5; do
    rt=${spec%%:*}
    knob=${spec##*:}
    nice ~/soil13/venv/bin/python samples_cv.py --build ../../../local/v14 --other ../../../local/v13 --read-types "$rt" \
        --knob "$knob" --variants "v14,+sample" --out ~/v14samples_$rt -t 6 > ~/v14samples_$rt.log 2>&1
done
echo finished > ~/v14samples_other.done
