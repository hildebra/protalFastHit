#!/bin/bash
# The strain stage as a unit that includes only what it uses (Options.h, Profiler.h and what they bring) instead of
# all of RunProtal.h: RunProtal.h's lines 2403-3734 (SharedSNP .. StrainWrapper2) behind the two aliases they need.
# Timed alone on one core, with and without the x86-64-v3 clones.
W=~/ctime
D=$W/split
H=$W/src/src/RunProtal.h
cd $D
{
    echo '#include <filesystem>'
    echo '#include <fstream>'
    echo '#include <iostream>'
    echo '#include <optional>'
    echo '#include <unistd.h>'
    echo '#include <sys/wait.h>'
    echo '#include "Options.h"'
    echo '#include "Profiler.h"'
    echo '#include "ProgressBar.h"'
    echo 'namespace protal {'
    echo '    using Profile = profiler::MicrobialProfile;'
    echo '    using Profiles = std::vector<Profile>;'
    sed -n '153,162p' $H     # MakeFolder
    sed -n '1872,1893p' $H   # MSAItemFor, MSAReleaseFor
    sed -n '2403,3734p' $H | sed -e 's/^    TaxidList ExtractTaxa/    static TaxidList ExtractTaxa/'
    echo '    namespace split {'
    echo '        void Strains(Options& options, Profiles& profiles, GenomeLoader& loader, taxonomy::IntTaxonomy& taxonomy,'
    echo '                     std::vector<uint32_t> msa_taxids, std::optional<profiler::TaxonFilterObj> filter) {'
    echo '            StrainWrapper2(options, profiles, loader, taxonomy, std::move(msa_taxids), std::move(filter));'
    echo '        }'
    echo '    }'
    echo '}'
} > lean_strains.cpp
BASE=$(cat base_cmd.txt)
taskset -c 3 $BASE -fsyntax-only lean_strains.cpp 2>&1 | grep -E 'error' | head -20
t0=$(date +%s.%N); taskset -c 3 nice -n 5 $BASE -fsyntax-only lean_strains.cpp 2>/dev/null; t1=$(date +%s.%N)
echo "lean_strains front end: $(echo "$t1 - $t0" | bc) s, preprocessed lines $($BASE -E lean_strains.cpp 2>/dev/null | wc -l)"
/usr/bin/time -f "lean_strains: %e s wall, %U s user" taskset -c 3 nice -n 5 $BASE -c lean_strains.cpp -o lean_strains.o
/usr/bin/time -f "lean_strains noclones: %e s wall, %U s user" taskset -c 3 nice -n 5 $BASE -DPROTAL_NO_CLONES -c lean_strains.cpp -o nc_lean_strains.o
for o in lean_strains.o nc_lean_strains.o; do
    echo "$o: text $(size $o | awk 'NR==2 {print $1}'), PrepareMAPQ $(nm -C $o | grep -c PrepareMAPQ), DecodeChunk $(nm -C $o | grep -c DecodeChunk), StrainWrapper2 $(nm -C $o | grep -c StrainWrapper2), GetMSAForTaxon $(nm -C $o | grep -c GetMSAForTaxon)"
done
