#!/bin/bash
# Splits ab6bf25's src/RunProtal.h into the stage units of src/Run/: every piece moved verbatim by line range, only the
# calls between the stages changed. Each range's first and last line are checked before it is cut.
set -euo pipefail
S=$(cd "$(dirname "$0")" && pwd)
W=C:/Users/hildebra/Documents/locDev/protal-split
P=$S/parts
O=$S/RunProtal.orig.h
git -C $W show ab6bf25:src/RunProtal.h > $O
D=$W/src/Run
mkdir -p $D

line() { sed -n "$1p" $O; }
expect() {  # line, text it must start with
    local got
    got=$(line $1)
    if [[ "$got" != "$2"* ]]; then echo "line $1: expected '$2', got '$got'" >&2; exit 1; fi
}
r() { sed -n "$1,$2p" $O; }

expect 54 '    static const std::string PROTAL_LOGO'; expect 72 '    }'
expect 160 '    // Alignments are written'; expect 168 '    }'
expect 170 '    // The index is built from --reference'; expect 242 '    }'
expect 244 "    // The threads a table's reading"; expect 353 '    }'
expect 355 "    // The database's models of kind"; expect 574 '    }'
expect 517 '    static build_settings::Settings DatabaseBuildSettings('
expect 576 "    // The species' NCBI taxids"; expect 824 '    }'
expect 826 '    // One line per sample'; expect 871 '    }'
expect 873 "    // Loads the database's index"; expect 904 '    }'
expect 906 '    // Whether a sample of the run is gated'; expect 913 '    }'
expect 915 '    // The seed expectation'; expect 951 '    }'
expect 953 '    // The grid level of the sparse index'; expect 1050 '    }'
expect 1052 "    // Whether a sample's reads come through"; expect 1063 '    }'
expect 1065 '    // --stage1_probe:'; expect 1107 '    }'
expect 1109 '    // Aligns the samples.'; expect 1650 '    }'
expect 1120 '        if (options.BuildMode()) {'; expect 1148 '        } else {'
expect 1655 '//     void ProfileWrapper2'; expect 1684 '//     }'
expect 1703 '    // The UnreportedTaxon of each sample'; expect 1737 '    }'
expect 1826 '    // Profiles the sample of slot idx'; expect 1869 '    }'
expect 1895 '    static void ProfileSample('; expect 2402 '    }'
expect 2404 '    static std::pair<bool,bool> SharedSNP'; expect 3733 '    }'
expect 3736 '    size_t GetTotalSystemMemory()'; expect 3741 '    }'
expect 3743 '    // Loads the PMML presence model'; expect 3785 '    }'
expect 3787 '    // Runs protal and returns'; expect 4311 '    }'; expect 4312 '}'
expect 3819 '        using AlignmentBenchmark = CoreBenchmark;'
expect 3841 '            // Each exits 8 on failure.'; expect 3880 '            return RunStatus::Get().Finish();'; expect 3881 '        }'
expect 3963 '        if (concurrent && run_alignment) {'; expect 3976 '        }'
expect 4225 '        // Untangle Template'; expect 4239 '        }'

# Load.cpp
{
    cat $P/Load.head.cpp
    { r 244 353; echo; r 576 824; echo; r 873 904; echo; r 953 1050; } | sed \
        -e 's/^    static void Load\([A-Za-z]*\)(Options const& options, GenomeLoader& genomes, std::ostream& out = std::cout) {$/    void Load\1(Options const\& options, GenomeLoader\& genomes, std::ostream\& out) {/' \
        -e 's/^    static std::optional<gene_neighbours::Table> ReadGeneNeighbours(Options const& options, std::ostream& out = std::cout) {$/    std::optional<gene_neighbours::Table> ReadGeneNeighbours(Options const\& options, std::ostream\& out) {/' \
        -e 's/^    static void BindGeneNeighbours(Options const& options, ProtalDB& db, gene_neighbours::Table table, std::ostream& out = std::cout) {$/    void BindGeneNeighbours(Options const\& options, ProtalDB\& db, gene_neighbours::Table table, std::ostream\& out) {/' \
        -e 's/^    static void RequireNcbiTaxids(/    void RequireNcbiTaxids(/' \
        -e 's/^    static bool GridWhole(/    bool GridWhole(/' \
        -e 's/^    static void CheckGridIndex(/    void CheckGridIndex(/' \
        -e 's/^    static std::unique_ptr<Seedmap> LoadIndex(/    std::unique_ptr<Seedmap> LoadIndex(/'
    cat $P/Load.tail.cpp
} > $D/Load.cpp

# Align.cpp
{
    cat $P/Align.head.cpp
    r 160 168; echo; r 826 871; echo; r 915 951; echo; r 1065 1107; echo
    { r 1109 1120; echo '            RunBuild(options, db);'; r 1148 1650; } | sed \
        -e 's/genera = profiler::GeneraOf(db.GetTaxonomy());/genera = TaxonomyGenera(db.GetTaxonomy());/'
    cat $P/Align.tail.cpp
} > $D/Align.cpp

# BuildStage.cpp
{
    cat $P/BuildStage.head.cpp
    r 170 242
    echo
    echo '    void RunBuild(Options& options, ProtalDB& db) {'
    echo '        const size_t mmer_size = 15;'
    echo '        const size_t kmer_size = 31;'
    echo
    r 1121 1146 | sed 's/^    //'
    echo '    }'
    echo
    echo '    int RunDbMaintenance(Options& options) {'
    r 3841 3880 | sed -e 's/^    //' \
        -e 's/^\( *\)LoadModel(db::DbFile::OnDisk(model.file), options.GetKnob(), &kind);/\1kind = ModelKindOfFile(db::DbFile::OnDisk(model.file), options.GetKnob());/'
    echo '    }'
    echo '}'
} > $D/BuildStage.cpp

# Profile.cpp
{
    cat $P/Profile.head.cpp
    r 1703 1737; echo; r 1826 1869; echo; r 1895 2402 | sed \
        -e 's/^    Profiles ProfileWrapper(ProfilingContext& ctx, ProfilingAhead\* ahead = nullptr) {$/    Profiles ProfileWrapper(ProfilingContext\& ctx, ProfilingAhead* ahead) {/'
    echo; r 1655 1684
    cat $P/Profile.tail.cpp
} > $D/Profile.cpp

# Strains.cpp
{
    cat $P/Strains.head.cpp
    r 2404 3733 | sed \
        -e '/^    using TaxidList = std::vector<uint32_t>;$/d' \
        -e 's/^    TaxidList ExtractTaxa(Profiles const& profiles, std::optional<profiler::TaxonFilterObj> filter= {}, size_t min_samples = 2) {$/    TaxidList ExtractTaxa(Profiles const\& profiles, std::optional<profiler::TaxonFilterObj> filter, size_t min_samples) {/' \
        -e 's/^    static std::vector<uint32_t> ResolveMSASpecies(/    std::vector<uint32_t> ResolveMSASpecies(/' \
        -e 's/^    static void StrainWrapper2(\(.*\), std::vector<uint32_t> msa_taxids = {}, std::optional<profiler::TaxonFilterObj> filter={}) {$/    void StrainWrapper2(\1, std::vector<uint32_t> msa_taxids, std::optional<profiler::TaxonFilterObj> filter) {/'
    echo '}'
} > $D/Strains.cpp

# Run.cpp
{
    cat $P/Run.head.cpp
    r 54 72; echo
    r 355 574 | sed -e 's/^    static build_settings::Settings DatabaseBuildSettings(/    build_settings::Settings DatabaseBuildSettings(/'
    echo
    r 3736 3741 | sed -e 's/^    size_t GetTotalSystemMemory()/    static size_t GetTotalSystemMemory()/'
    echo
    r 3743 3785
    echo
    echo '    std::string ModelKindOfFile(db::DbFile const& file, double knob) {'
    echo '        std::string kind;'
    echo '        LoadModel(file, knob, &kind);'
    echo '        return kind;'
    echo '    }'
    echo
    { r 3787 3818; r 3820 3840
      echo '            return RunDbMaintenance(options);  // each exits 8 on failure'
      r 3881 3962
      echo '        if (concurrent && run_alignment) StartIndexLoading(options, db);'
      r 3977 4224
      echo '        if (run_alignment) RunAlignment(options, db, sample_done);'
      r 4240 4312; } | sed \
        -e 's/^    static int Run(int argc, char \*argv\[\]) {$/    int Run(int argc, char *argv[]) {/' \
        -e 's/^        std::unique_ptr<ProfilingAhead> ahead;$/        ProfilingAheadPtr ahead;/' \
        -e 's/ahead = std::make_unique<ProfilingAhead>(/ahead = MakeProfilingAhead(/' \
        -e 's/ahead->Submit(/SubmitProfilingAhead(*ahead, /'
} > $D/Run.cpp

# The checks that each edit took
check() {  # file, pattern, expected count
    local n
    n=$(grep -c -- "$2" "$1" || true)
    if [[ "$n" != "$3" ]]; then echo "$(basename $1): '$2' found $n times, expected $3" >&2; exit 1; fi
}
check $D/Load.cpp 'std::ostream& out = std::cout' 0
check $D/Load.cpp '^    static' 6
check $D/Align.cpp 'RunBuild(options, db);' 1
check $D/Align.cpp 'TaxonomyGenera(' 1
check $D/Align.cpp 'profiler::' 0
check $D/BuildStage.cpp 'ModelKindOfFile(' 1
check $D/BuildStage.cpp 'LoadModel' 0
check $D/Profile.cpp 'ProfilingAhead\* ahead = nullptr' 0
check $D/Strains.cpp 'using TaxidList' 0
check $D/Strains.cpp '^    TaxidList ExtractTaxa(Profiles const& profiles, std::optional<profiler::TaxonFilterObj> filter, size_t min_samples) {' 1
check $D/Strains.cpp '^    std::vector<uint32_t> ResolveMSASpecies(' 1
check $D/Strains.cpp '^    void StrainWrapper2(' 1
check $D/Run.cpp '^    int Run(int argc, char \*argv\[\]) {' 1
check $D/Run.cpp 'ProfilingAheadPtr ahead;' 1
check $D/Run.cpp 'MakeProfilingAhead(' 1
check $D/Run.cpp 'SubmitProfilingAhead(\*ahead, ' 1
check $D/Run.cpp 'RunWrapper' 0
check $D/Run.cpp 'protal::build::[A-Z]' 2
check $D/Run.cpp 'AlignmentBenchmark' 0
check $D/Run.cpp '^    build_settings::Settings DatabaseBuildSettings(' 1
wc -l $D/*.cpp $D/*.h
