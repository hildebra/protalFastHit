#!/bin/bash
# Prototype: protal::Run's stages in translation units of their own, every unit including all of RunProtal.h
# (no header changes). Each unit's -O3 compile timed, the units compiled 4 at a time, then linked as protal.
set -e
S=$(cd "$(dirname "$0")" && pwd)  # SplitDecls.h is next to this script
W=~/ctime
B=$W/build
H=$W/src/src/RunProtal.h
D=$W/split
rm -rf $D
mkdir -p $D
cp $S/SplitDecls.h $D/
# Bottom to top; every address is a line of the original file.
sed -e '4298s/StrainWrapper2(/split::Strains(/' \
    -e '4247s/ProfileWrapper(/split::Profile(/' \
    -e '4225,4238c\        if (run_alignment) split::Align(options, db, sample_done);' \
    -e '4214s/std::make_unique<ProfilingAhead>(/split::MakeAhead(/' \
    -e '3840,3880d' \
    -e '3839s/) {$/) return split::BuildMaint(options);/' \
    -e '3735i #include "SplitDecls.h"' \
    -e '3735s/^    size_t/    static size_t/' \
    -e '2625s/^    TaxidList/    static TaxidList/' \
    -e '2398s/^    Profiles/    static Profiles/' \
    -e '2338s/^    Profiles/    static Profiles/' \
    -e '1120,1145d' \
    -e '1119a\            split::BuildIndex(options, db);' \
    -e '3803,3836s/protal::Options::OptionsFromArguments(argc, argv)/split::ParseOptions(argc, argv)/' \
    -e '3803,3836s/options\.PrintHelp(true);/split::Help(options, true);/' \
    -e '3803,3836s/options\.PrintHelp();/split::Help(options, false);/' \
    -e '3803,3836s/options\.ToString()/split::OptionsText(options)/' \
    -e '1111i\    namespace split { void BuildIndex(Options& options, ProtalDB& db); }' \
    $H > $D/RunProtal.h
cp $W/src/src/main.cpp $D/main.cpp
{
    echo '#include "RunProtal.h"'
    echo 'namespace protal::split {'
    echo '    void Align(Options& options, ProtalDB& db, std::function<void(size_t)> const& sample_done) {'
    echo '        using AlignmentBenchmark = CoreBenchmark;'
    sed -n '4226,4237p' $H
    echo '    }'
    echo '}'
} > $D/split_align.cpp
{
    echo '#include "RunProtal.h"'
    echo 'namespace protal::split {'
    echo '    Profiles Profile(ProfilingContext& ctx, ProfilingAhead* ahead) { return ProfileWrapper(ctx, ahead); }'
    echo '    std::unique_ptr<ProfilingAhead> MakeAhead(ProfilingContext& ctx, size_t threads) {'
    echo '        return std::make_unique<ProfilingAhead>(ctx, threads);'
    echo '    }'
    echo '}'
} > $D/split_profile.cpp
{
    echo '#include "RunProtal.h"'
    echo 'namespace protal::split {'
    echo '    void Strains(Options& options, Profiles& profiles, GenomeLoader& loader, taxonomy::IntTaxonomy& taxonomy,'
    echo '                 std::vector<uint32_t> msa_taxids, std::optional<profiler::TaxonFilterObj> filter) {'
    echo '        StrainWrapper2(options, profiles, loader, taxonomy, std::move(msa_taxids), std::move(filter));'
    echo '    }'
    echo '}'
} > $D/split_strains.cpp
{
    echo '#include "RunProtal.h"'
    echo 'namespace protal::split {'
    echo '    int BuildMaint(Options& options) {'
    sed -n '3840,3879p' $H
    echo '    }'
    echo '}'
} > $D/split_buildmaint.cpp
{
    echo '#include "RunProtal.h"'
    echo 'namespace protal::split {'
    echo '    void BuildIndex(Options& options, ProtalDB& db) {'
    echo '        const size_t mmer_size = 15;'
    echo '        const size_t kmer_size = 31;'
    sed -n '1120,1145p' $H
    echo '    }'
    echo '}'
} > $D/split_buildindex.cpp
{
    echo '#include "RunProtal.h"'
    echo 'namespace protal::split {'
    echo '    Options ParseOptions(int argc, char* argv[]) { return protal::Options::OptionsFromArguments(argc, argv); }'
    echo '    void Help(Options& options, bool dev) { options.PrintHelp(dev); }'
    echo '    std::string OptionsText(Options const& options) { return options.ToString(); }'
    echo '}'
} > $D/split_options.cpp
diff $H $D/RunProtal.h | head -80 || true
CMD=$(grep -o '"command": "[^"]*protal.dir/src/main.cpp[^"]*"' $B/compile_commands.json | head -1 | sed 's/^"command": "//; s/"$//')
BASE=${CMD% -o *}
echo "$BASE" > $D/base_cmd.txt
cd $D
for u in main split_align split_buildindex split_options split_profile split_strains split_buildmaint; do
    taskset -c 0-3 $BASE -fsyntax-only $u.cpp
done
echo "syntax ok"
uptime
compile() {
    local u=$1
    /usr/bin/time -f "$u: %e s wall, %U s user, %M KB" taskset -c 0-3 nice -n 5 $(cat base_cmd.txt) -ftime-report -c $u.cpp -o $u.o 2> ft_$u.txt
    tail -1 ft_$u.txt
}
export -f compile
start=$(date +%s.%N)
printf '%s\n' main split_align split_buildindex split_options split_profile split_strains split_buildmaint | xargs -P4 -I{} bash -c 'compile {}'
end=$(date +%s.%N)
echo "all units, 4 at a time: $(echo "$end - $start" | bc) s wall"
uptime
for u in main split_align split_buildindex split_options split_profile split_strains split_buildmaint; do
    echo "$u $(grep -E '^ TOTAL' ft_$u.txt) text=$(size $u.o | awk 'NR==2 {print $1}')"
done
# Link as protal: the build's own link command, main.cpp.o replaced by the five units.
cd $B
LINK=$(ninja -t commands protal | tail -1)
echo "$LINK" > $D/link_cmd.txt
LINK=${LINK/CMakeFiles\/protal.dir\/src\/main.cpp.o/$D/main.o $D/split_align.o $D/split_profile.o $D/split_strains.o $D/split_buildmaint.o $D/split_buildindex.o $D/split_options.o}
LINK=${LINK/-o protal /-o $D/protal_split }
/usr/bin/time -f "link: %e s wall" bash -c "$LINK" && ls -la $D/protal_split && $D/protal_split --version | head -3
