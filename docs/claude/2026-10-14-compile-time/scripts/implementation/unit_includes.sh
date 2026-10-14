#!/bin/bash
# Which of protal's big headers each unit of the split build includes (g++ -H on the unit's own compile command), and
# its preprocessed size.
B=~/split-build/build
cd $B
for u in src/main.cpp src/Run/Run.cpp src/Run/Load.cpp src/Run/Align.cpp src/Run/BuildStage.cpp src/Run/Profile.cpp src/Run/Strains.cpp src/Options.cpp; do
    cmd=$(grep -o "\"command\": \"[^\"]*/$u\"" compile_commands.json | head -1 | sed 's/^"command": "//; s/"$//')
    cmd=${cmd%% -o *}
    dir=$(grep -B1 "\"command\": \"[^\"]*/$u\"" compile_commands.json | grep -o '"directory": "[^"]*"' | head -1 | sed 's/"directory": "//; s/"$//')
    lines=$(cd ${dir:-$B} && $cmd -E ~/split-build/src/$u 2>/dev/null | wc -l)
    big=$(cd ${dir:-$B} && $cmd -H -E ~/split-build/src/$u -o /dev/null 2>&1 | grep -oE '/(Options|GenomeLoader|Profiler|Classify|Build|Seedmap|IndexCodec|KmerLookup|AlignmentStrategy|ChainAnchorFinder|Strain|cxxopts|cPMML|ThreadedGzStream)\.(h|hpp)$' | sort -u | sed 's#^/##' | tr '\n' ' ')
    printf '%-24s %7d lines  %s\n' $u $lines "$big"
done
