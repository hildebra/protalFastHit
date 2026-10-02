// Unit tests for the database and input files protal reads: malformed files must stop with a clear
// message (exit 8) instead of being half-read, and the index must know which reference it was built for.
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h>
#include "Hash/Seedmap.h"
#include "Profiling/Profiler.h"
#include "SequenceUtils/GenomeLoader.h"
#include "SequenceUtils/KmerUtils.h"
#include "Taxonomy/Taxonomy.h"
#include "Utilities/ReferenceFingerprint.h"
#include "SequenceUtils/ReadTypeDetection.h"
#include "gzstream/gzstream.h"

namespace fs = std::filesystem;

namespace {
    struct ScratchDir {
        fs::path path;
        ScratchDir() {
            path = fs::temp_directory_path() / ("protal input test " + std::to_string(::getpid()));
            fs::create_directories(path);
        }
        ~ScratchDir() { fs::remove_all(path); }

        std::string Write(std::string const& name, std::string const& content) const {
            auto file = path / name;
            std::ofstream(file, std::ios::binary) << content;
            return file.string();
        }
    };

    // Two genera with two and one species, as the mini DB lays them out.
    constexpr char kTaxonomy[] =
            "id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n"
            "1\t6\t0\ts__Mockella alpha\tspecies\t7\tGCF_1\n"
            "2\t6\t0\ts__Mockella beta\tspecies\t7\tGCF_2\n"
            "3\t7\t0\ts__Fakibacter gamma\tspecies\t7\tGCF_3\n"
            "4\t4\t0\troot\tno rank\t0\t\n"
            "5\t4\t0\td__Bacteria\tdomain\t1\t\n"
            "6\t5\t0\tg__Mockella\tgenus\t6\t\n"
            "7\t5\t0\tg__Fakibacter\tgenus\t6\t\n";

    void LoadTaxonomy(std::string const& path) { protal::taxonomy::IntTaxonomy taxonomy(path); }
}

TEST(Taxonomy, LoadsWithOrWithoutHeaderAndWithCRLF) {
    ScratchDir dir;
    protal::taxonomy::IntTaxonomy with_header(dir.Write("a.dmp", kTaxonomy));
    EXPECT_EQ(with_header.Get("s__Mockella beta"), 2u);
    EXPECT_EQ(with_header.root_id, 4u);

    std::string body = kTaxonomy;
    body = body.substr(body.find('\n') + 1);
    protal::taxonomy::IntTaxonomy without_header(dir.Write("b.dmp", body));
    EXPECT_EQ(without_header.map.size(), with_header.map.size());

    std::string crlf;
    for (char c : std::string(kTaxonomy)) crlf += c == '\n' ? std::string("\r\n") : std::string(1, c);
    protal::taxonomy::IntTaxonomy windows(dir.Write("c.dmp", crlf + "\r\n"));
    EXPECT_EQ(windows.map.at(3).rank, "species");
    EXPECT_EQ(windows.map.at(3).rep_genome, "GCF_3");
}

TEST(Taxonomy, RejectsMalformedFiles) {
    ScratchDir dir;
    std::string header = "id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n";
    EXPECT_EXIT(LoadTaxonomy((dir.path / "none.dmp").string()),
                testing::ExitedWithCode(8), "cannot open");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("empty.dmp", header)), testing::ExitedWithCode(8), "no taxa");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("short.dmp", header + "1\t1\t0\troot\tno rank\n")),
                testing::ExitedWithCode(8), "line 2: expected at least 6");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("text.dmp", header + "1\tone\t0\troot\tno rank\t0\n")),
                testing::ExitedWithCode(8), "must be integers");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("twice.dmp", std::string(kTaxonomy) + "2\t6\t0\ts__Other\tspecies\t7\t\n")),
                testing::ExitedWithCode(8), "id 2 is defined twice");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("orphan.dmp", std::string(kTaxonomy) + "8\t9\t0\ts__Lost\tspecies\t7\t\n")),
                testing::ExitedWithCode(8), "parent id 9 is used but never defined");
    // Lineages are walked up to the root, the taxon that is its own parent: a cycle, or a root
    // given a parent, would make the walk go on for ever.
    EXPECT_EXIT(LoadTaxonomy(dir.Write("cycle.dmp", std::string(kTaxonomy) + "8\t9\t0\ts__Round\tspecies\t7\t\n9\t8\t0\tg__Round\tgenus\t6\t\n")),
                testing::ExitedWithCode(8), "never reaches the root .a cycle through id [89]");
    std::string rooted_below = kTaxonomy;
    rooted_below.replace(rooted_below.find("4\t4\t0\troot"), 3, "4\t5");
    EXPECT_EXIT(LoadTaxonomy(dir.Write("rootless.dmp", rooted_below)), testing::ExitedWithCode(8), "never reaches the root");
    // Truth files and --msa_species name taxa, so a name must be one taxon's.
    EXPECT_EXIT(LoadTaxonomy(dir.Write("same.dmp", std::string(kTaxonomy) + "8\t6\t0\ts__Mockella alpha\tspecies\t7\t\n")),
                testing::ExitedWithCode(8), "line 9: the name 's__Mockella alpha' is given to ids 1 and 8; names must be unique");
}

TEST(Truth, ReadsLineagesAndTaxids) {
    ScratchDir dir;
    protal::taxonomy::IntTaxonomy taxonomy(dir.Write("t.dmp", kTaxonomy));

    // simulate_metagenomes writes the lineage in a later column; any field may hold it.
    auto lineage = protal::GetTruth(dir.Write("lineage.tsv",
            "genome\tlineage\treads\n"
            "GCF_1\td__Bacteria;g__Mockella;s__Mockella alpha\t100\r\n"
            "\n"
            "# a comment\n"
            "GCF_3\td__Bacteria;g__Fakibacter;s__Fakibacter gamma\t50\n"), taxonomy);
    EXPECT_EQ(lineage.size(), 2u);
    EXPECT_TRUE(lineage.contains(1));
    EXPECT_TRUE(lineage.contains(3));

    // Single-column species names and internal taxids.
    auto names = protal::GetTruth(dir.Write("names.tsv", "s__Mockella beta\n2\n3\t0.4\n"), taxonomy);
    EXPECT_EQ(names.size(), 2u);
    EXPECT_TRUE(names.contains(2));
    EXPECT_TRUE(names.contains(3));
}

TEST(Truth, SkipsLinesWithoutAKnownSpecies) {
    ScratchDir dir;
    protal::taxonomy::IntTaxonomy taxonomy(dir.Write("t.dmp", kTaxonomy));
    testing::internal::CaptureStderr();
    auto truth = protal::GetTruth(dir.Write("partial.tsv",
            "GCF_9\td__Bacteria;g__Mockella;s__\n"           // unclassified at species level
            "GCF_8\td__Bacteria;s__Unknown species\n"
            "99\n"                                           // taxid not in the taxonomy
            "just text\n"
            "GCF_2\td__Bacteria;g__Mockella;s__Mockella beta\n"), taxonomy);
    auto log = testing::internal::GetCapturedStderr();
    EXPECT_EQ(truth.size(), 1u);
    EXPECT_TRUE(truth.contains(2));
    EXPECT_NE(log.find("No species classification"), std::string::npos);
    EXPECT_NE(log.find("Species unknown to taxonomy: s__Unknown species"), std::string::npos);
    EXPECT_NE(log.find("taxid 99 is not in the taxonomy"), std::string::npos);
    EXPECT_NE(log.find("line 4: no lineage or taxid"), std::string::npos);
}

TEST(Truth, MissingFileStops) {
    ScratchDir dir;
    protal::taxonomy::IntTaxonomy taxonomy(dir.Write("t.dmp", kTaxonomy));
    EXPECT_EXIT(protal::GetTruth((dir.path / "none.tsv").string(), taxonomy), testing::ExitedWithCode(8), "Cannot read truth file");
}

namespace {
    // Two genes of taxon 1 and one of taxon 2; the map holds the byte range of each sequence line.
    struct Reference {
        std::string fna;
        std::string map;
        Reference() {
            std::ostringstream fna_os, map_os;
            std::pair<int, int> genes[] = { { 1, 1 }, { 1, 2 }, { 2, 1 } };
            std::string sequences[] = { "ACGTACGTAA", "ccggttaa", "GGGGCCCCAT" };
            for (int i = 0; i < 3; i++) {
                fna_os << '>' << genes[i].first << '_' << genes[i].second << '\n';
                size_t start = fna_os.str().size();
                fna_os << sequences[i] << '\n';
                map_os << genes[i].first << '\t' << genes[i].second << '\t' << start << '\t' << start + sequences[i].size() << '\n';
            }
            fna = fna_os.str();
            map = map_os.str();
        }
    };

    void LoadMap(std::string const& fna, std::string const& map) { protal::GenomeLoader loader(fna, map); }
}

TEST(ReferenceMap, LoadsGenesAndUppercasesThem) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    protal::GenomeLoader loader(fna, dir.Write("reference.map", ref.map + "\n"));
    EXPECT_EQ(loader.GeneCount(), 3u);
    EXPECT_TRUE(loader.HasGene(1, 2));
    EXPECT_FALSE(loader.HasGene(2, 2));
    EXPECT_FALSE(loader.HasGene(3, 1));
    EXPECT_EQ(loader.GeneLength(1, 2), 8u);

    auto& gene = loader.GetGenome(1).GetGeneOMP(2);
    EXPECT_EQ(gene.Sequence(), "CCGGTTAA");
}

// LoadAllGenomes puts the genes not yet loaded into one arena; a genome loaded gene by gene
// before keeps its sequences.
TEST(ReferenceMap, LoadAllGenomesAfterAGenomeWasLoadedOnItsOwn) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    for (int threads : { 1, 3 }) {
        protal::GenomeLoader loader(fna, map);
        EXPECT_EQ(loader.GetGenome(1).GetGeneOMP(2).Sequence(), "CCGGTTAA");
        loader.LoadAllGenomes(threads);
        EXPECT_TRUE(loader.AllGenomesLoaded());
        EXPECT_EQ(loader.GetGenome(1).GetGene(1).Sequence(), "ACGTACGTAA") << threads << " threads";
        EXPECT_EQ(loader.GetGenome(1).GetGene(2).Sequence(), "CCGGTTAA") << threads << " threads";
        EXPECT_EQ(loader.GetGenome(2).GetGene(1).Sequence(), "GGGGCCCCAT") << threads << " threads";

        protal::GenomeLoader fresh(fna, map);
        fresh.LoadAllGenomes(threads);
        EXPECT_EQ(fresh.GetGenome(1).GetGene(2).Sequence(), "CCGGTTAA") << threads << " threads";
        EXPECT_EQ(fresh.GetGenome(2).GetGene(1).Sequence(), "GGGGCCCCAT") << threads << " threads";
    }
}

TEST(ReferenceMap, RejectsMalformedLines) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto size = std::to_string(ref.fna.size());
    EXPECT_EXIT(LoadMap(fna, dir.Write("m1", ref.map + "3\t1\t5\n")), testing::ExitedWithCode(8), "line 4: expected 4 tab-separated columns, found 3");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m2", ref.map + "3\t1\t5\tx\n")), testing::ExitedWithCode(8), "column 4 is not");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m3", ref.map + "0\t1\t5\t8\n")), testing::ExitedWithCode(8), "taxid must be between 1 and");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m4", ref.map + "3\t0\t5\t8\n")), testing::ExitedWithCode(8), "gene id must be between 1 and");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m5", ref.map + "3\t1\t8\t8\n")), testing::ExitedWithCode(8), "end byte must be after start byte");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m6", ref.map + "3\t1\t5\t" + size + "0\n")), testing::ExitedWithCode(8), "is past the end of");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m7", ref.map + "1\t2\t5\t8\n")), testing::ExitedWithCode(8), "gene 1_2 is listed twice");
    EXPECT_EXIT(LoadMap(fna, (dir.path / "none.map").string()), testing::ExitedWithCode(8), "cannot open the file");
    EXPECT_EXIT(LoadMap(fna, dir.Write("m8", "\n")), testing::ExitedWithCode(8), "the file lists no genes");
}

TEST(ReferenceFingerprint, TracksMapContentAndReferenceSize) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    auto original = protal::ReferenceFingerprint::Of(map, fna);
    EXPECT_EQ(original.fna_size, ref.fna.size());
    EXPECT_EQ(original, protal::ReferenceFingerprint::Of(map, fna));

    dir.Write("reference.map", ref.map + "3\t1\t5\t8\n");
    EXPECT_NE(original, protal::ReferenceFingerprint::Of(map, fna));
    dir.Write("reference.map", ref.map);
    dir.Write("reference.fna", ref.fna + "A");
    EXPECT_NE(original, protal::ReferenceFingerprint::Of(map, fna));
}

TEST(IndexHeader, ReferenceFingerprintRoundTrips) {
    protal::ReferenceFingerprint fingerprint{ 0x0123456789abcdefULL, 386271 };
    std::stringstream with;
    {
        protal::Seedmap built;
        EXPECT_FALSE(built.HasReferenceFingerprint());
        built.SetReferenceFingerprint(fingerprint);
        built.SaveHeader(with);
    }
    protal::Seedmap loaded;
    loaded.LoadHeader(with);
    ASSERT_TRUE(loaded.HasReferenceFingerprint());
    EXPECT_EQ(loaded.GetReferenceFingerprint(), fingerprint);
    EXPECT_TRUE(loaded.UsesFullSyncmerMask());
    EXPECT_NE(loaded.FeatureDescription().find(", reference fingerprint"), std::string::npos);

    // Indices written before the fingerprint existed still load, without one.
    std::stringstream without;
    {
        protal::Seedmap built;
        built.SaveHeader(without);
    }
    protal::Seedmap old;
    old.LoadHeader(without);
    EXPECT_FALSE(old.HasReferenceFingerprint());
    EXPECT_NE(old.FeatureDescription().find("no reference fingerprint"), std::string::npos);
}

TEST(KmerUtils, LowercaseBasesEncodeLikeUppercase) {
    for (auto [upper, lower] : { std::pair{ 'A', 'a' }, { 'C', 'c' }, { 'G', 'g' }, { 'T', 't' } }) {
        EXPECT_EQ(KmerUtils::BaseToInt(lower), KmerUtils::BaseToInt(upper));
        EXPECT_EQ(KmerUtils::BaseToIntC(lower), KmerUtils::BaseToIntC(upper));
        EXPECT_EQ(KmerUtils::BaseToInt(lower, 0), KmerUtils::BaseToInt(upper, 0));
        EXPECT_EQ(KmerUtils::BaseToIntC(lower, 0), KmerUtils::BaseToIntC(upper, 0));
    }
    EXPECT_EQ(KmerUtils::BaseToInt('n'), 4u);
}

TEST(UniqueKmers, LoadsCountsAndRejectsMalformedLines) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    {
        protal::GenomeLoader loader(fna, map);
        loader.LoadUniqueKmers(dir.Write("u.tsv", "1\t1\t3\t0.3\t0\t0\t0\t0\t10\r\n\n1\t2\t0\t0\t0\t0\t0\t0\t8\n"));
        EXPECT_TRUE(loader.GetGenome(1).IsGeneHittable(1));
        EXPECT_FALSE(loader.GetGenome(1).IsGeneHittable(2));
    }
    auto load = [&](std::string const& content) {
        protal::GenomeLoader loader(fna, map);
        loader.LoadUniqueKmers(dir.Write("bad.tsv", content));
    };
    EXPECT_EXIT(load("1\t1\t3\n"), testing::ExitedWithCode(8), "line 1: expected 9 tab-separated columns, found 3");
    EXPECT_EXIT(load("1\t1\tx\t0\t0\t0\t0\t0\t10\n"), testing::ExitedWithCode(8), "column 3 is not");
    EXPECT_EXIT(load("7\t1\t3\t0.3\t0\t0\t0\t0\t10\n"), testing::ExitedWithCode(8), "gene 7_1 is not in reference.map");
    EXPECT_EXIT(load(""), testing::ExitedWithCode(8), "the file lists no genes .rebuild the database with --build.");
}

namespace {
    // A reference.map of 600 taxa x 100 genes (~1.7 MB: several chunks of the parallel parser) over a
    // sparse reference.fna, and a unique_kmers.tsv for it. `edit` changes lines (0-based) by index.
    struct BigTables {
        std::string fna, map, unique;
        size_t genes = 60000;
        BigTables(ScratchDir const& dir, std::map<size_t, std::string> const& map_edits = {},
                  std::map<size_t, std::string> const& unique_edits = {}) {
            std::string m, u;
            uint64_t position = 0;
            size_t line = 0;
            for (int taxid = 1; taxid <= 600; taxid++) {
                for (int gene = 1; gene <= 100; gene++, line++) {
                    uint64_t const length = 900 + (taxid * 7 + gene) % 300;
                    std::string map_line = std::to_string(taxid) + '\t' + std::to_string(gene) + '\t' + std::to_string(position) + '\t' + std::to_string(position + length);
                    std::string unique_line = std::to_string(taxid) + '\t' + std::to_string(gene) + '\t' + std::to_string(gene % 3) + "\t0.1\t" +
                                              std::to_string(gene % 2) + "\t0.1\t1\t0.1\t" + std::to_string(length - 30);
                    if (auto it = map_edits.find(line); it != map_edits.end()) map_line = it->second;
                    if (auto it = unique_edits.find(line); it != unique_edits.end()) unique_line = it->second;
                    m += map_line + '\n';
                    u += unique_line + '\n';
                    position += length + 1;
                }
            }
            fna = dir.Write("big.fna", "");
            std::filesystem::resize_file(fna, position);
            map = dir.Write("big.map", m);
            unique = dir.Write("big_unique.tsv", u);
        }
    };

    protal::GenomeLoader LoadBig(BigTables const& t, int threads) {
        protal::GenomeLoader loader(protal::db::DbFile::OnDisk(t.fna), protal::db::DbFile::OnDisk(t.map), threads);
        loader.LoadUniqueKmers(t.unique, threads);
        return loader;
    }
}

TEST(GeneTables, ParsedInParallelAsLineByLine) {
    ScratchDir dir;
    BigTables tables(dir);
    ASSERT_GT(std::filesystem::file_size(tables.map), size_t{1} << 20);  // several chunks
    auto one = LoadBig(tables, 1);
    auto eight = LoadBig(tables, 8);
    EXPECT_EQ(one.GeneCount(), tables.genes);
    EXPECT_EQ(eight.GeneCount(), tables.genes);
    for (int taxid : { 1, 299, 600 }) {
        for (int gene : { 1, 50, 100 }) {
            EXPECT_EQ(eight.GeneLength(taxid, gene), one.GeneLength(taxid, gene));
            EXPECT_EQ(eight.GetGenome(taxid).GetGene(gene).GetStartByte(), one.GetGenome(taxid).GetGene(gene).GetStartByte());
            EXPECT_EQ(eight.GetGenome(taxid).GetGene(gene).GetUniqueKmerCounts(), one.GetGenome(taxid).GetGene(gene).GetUniqueKmerCounts());
            EXPECT_EQ(eight.GetGenome(taxid).IsGeneHittable(gene), one.GetGenome(taxid).IsGeneHittable(gene));
        }
        EXPECT_EQ(eight.GetGenome(taxid).GetUniqueKmerCounts(), one.GetGenome(taxid).GetUniqueKmerCounts());
    }
}

TEST(GeneTables, TheFirstProblemInTheFileIsReportedWithItsLine) {
    ScratchDir dir;
    // A problem deep in the file (a late chunk) keeps its line number.
    BigTables late(dir, { { 55000, "7\t8\t9" } });
    EXPECT_EXIT(LoadBig(late, 8), testing::ExitedWithCode(8), "line 55001: expected 4 tab-separated columns, found 3");
    // A gene listed twice early wins over a malformed line later, and the other way round.
    BigTables duplicate_first(dir, { { 1000, "1\t1\t0\t10" }, { 50000, "x" } });
    EXPECT_EXIT(LoadBig(duplicate_first, 8), testing::ExitedWithCode(8), "line 1001: gene 1_1 is listed twice");
    BigTables malformed_first(dir, { { 1000, "1\t1\tx\t10" }, { 50000, "1\t1\t0\t10" } });
    EXPECT_EXIT(LoadBig(malformed_first, 8), testing::ExitedWithCode(8), "line 1001: column 3 is not a non-negative integer");
    BigTables huge(dir, { { 3, "1\t4\t0\t99999999999999999999999" } });
    EXPECT_EXIT(LoadBig(huge, 8), testing::ExitedWithCode(8), "line 4: column 4 is too large");
    // unique_kmers.tsv likewise.
    BigTables unknown_gene(dir, {}, { { 40000, "1\t700\t1\t0.1\t0\t0.1\t1\t0.1\t9" } });
    EXPECT_EXIT(LoadBig(unknown_gene, 8), testing::ExitedWithCode(8), "line 40001: gene 1_700 is not in reference.map");
    BigTables short_line(dir, {}, { { 59999, "1\t2" } });
    EXPECT_EXIT(LoadBig(short_line, 8), testing::ExitedWithCode(8), "line 60000: expected 9 tab-separated columns, found 2");
}

TEST(ModelFeatures, NormalizedFeaturesOfATaxon) {
    // Genome 1 has two hittable genes of 10 and 8 bp. Three mates of two fragments land on gene 1_1,
    // one with an insertion.
    ScratchDir dir;
    Reference ref;
    protal::GenomeLoader loader(dir.Write("reference.fna", ref.fna), dir.Write("reference.map", ref.map));
    loader.LoadAllGenomes();
    protal::profiler::MicrobialProfile profile(loader);
    auto sam = [](std::string seq, std::string cigar) {
        protal::SamEntry s;
        s.m_qname = "r";
        s.m_rname = "1_1";
        s.m_pos = 1;
        s.m_mapq = 60;
        s.m_qual = std::string(seq.size(), 'I');
        s.m_seq = std::move(seq);
        s.m_cigar = std::move(cigar);
        return s;
    };
    auto exact = sam("ACGTACGTAA", "10M");
    auto insertion = sam("ACGTAGCGTAA", "5M1I5M");  // identity 10/11
    ASSERT_TRUE(profile.AddSam(1, 1, exact, 1.0, true, 0));
    ASSERT_TRUE(profile.AddSam(1, 1, exact, 1.0, true, 0));  // its mate: the same fragment
    ASSERT_TRUE(profile.AddSam(1, 1, insertion, 1.0, true, 1));
    auto const& taxon = profile.GetTaxa().at(1);

    EXPECT_EQ(taxon.TotalHits(), 3u);
    EXPECT_EQ(taxon.Fragments(), 2u);
    EXPECT_NEAR(taxon.BaseIdentity(), (10 + 10 + 10 * 10.0 / 11) / 30, 1e-9);
    EXPECT_NEAR(taxon.TopIdentity(), 1.0, 1e-6);
    EXPECT_NEAR(taxon.HitGeneFraction(), 0.5, 1e-12);
    // Two fragments over genes of 10 and 8 bp: 1 - (8/18)^2 + 1 - (10/18)^2 genes expected, 1 hit.
    double const expected_genes = 2 - std::pow(8.0 / 18, 2) - std::pow(10.0 / 18, 2);
    EXPECT_NEAR(taxon.GenePresenceRatio(), 1 / expected_genes, 1e-9);
    // Fragments per gene 2 and 0 against 2 * 10/18 and 2 * 8/18.
    double const e1 = 2 * 10.0 / 18, e2 = 2 * 8.0 / 18;
    EXPECT_NEAR(taxon.GeneDispersion(), (2 - e1) * (2 - e1) / e1 + e2, 1e-9);
}

TEST(UniqueKmers, AGenomeWithoutUniqueKmersHasNoHittableGene) {
    ScratchDir dir;
    Reference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);

    protal::GenomeLoader without(fna, map);  // no unique_kmers.tsv: every gene counts
    EXPECT_TRUE(without.GetGenome(2).IsGeneHittable(1));
    EXPECT_EQ(without.GetGenome(2).GeneNum(), 1u);

    // Genome 2's gene has no unique k-mers, and genome 1 has one with and one without.
    protal::GenomeLoader loader(fna, map);
    loader.LoadUniqueKmers(dir.Write("u.tsv", "1\t1\t3\t0.3\t0\t0\t0\t0\t10\n1\t2\t0\t0\t0\t0\t0\t0\t8\n2\t1\t0\t0\t0\t0\t0\t0\t8\n"));
    EXPECT_EQ(loader.GetGenome(1).GetHittableGenes(), std::vector<uint32_t>{ 1 });
    EXPECT_FALSE(loader.GetGenome(2).IsGeneHittable(1));
    EXPECT_TRUE(loader.GetGenome(2).GetHittableGenes().empty());
    EXPECT_EQ(loader.GetGenome(2).GeneNum(), 0u);
}

namespace {
    // n reads of `length` random bases in FASTQ (or FASTA without quality), named name(i), every base of quality q
    // (Phred+33 character) but every tenth of quality low where low is given.
    std::string Reads(size_t n, size_t length, char q, std::function<std::string(size_t)> const& name, bool fasta = false, char low = 0) {
        std::mt19937 rng(5);
        std::string out;
        for (size_t i = 0; i < n; i++) {
            std::string seq(length, 'A'), qual(length, q);
            for (auto& c : seq) c = "ACGT"[rng() % 4];
            if (low) for (size_t k = 0; k < length; k += 10) qual[k] = low;
            out += (fasta ? ">" : "@") + name(i) + "\n" + seq + "\n";
            if (!fasta) out += "+\n" + qual + "\n";
        }
        return out;
    }
    std::string Plain(size_t i) { return "read" + std::to_string(i); }

    std::optional<protal::ReadTypeGuess> Guess(std::string const& path) {
        return protal::GuessReadType(protal::SampleReads(path), 1000);
    }
}

// The kind of reads in a single read file, from its first reads (ReadTypeDetection.h).
TEST(ReadTypeDetection, ShortAndLongReadsByLength) {
    ScratchDir dir;
    auto const short_reads = Guess(dir.Write("short.fq", Reads(300, 150, 'I', Plain)));
    ASSERT_TRUE(short_reads);
    EXPECT_EQ(short_reads->type, protal::ReadType::Single);
    // A stray long read among short ones leaves them short reads (the reader then stops at it).
    auto const stray = Guess(dir.Write("stray.fq", Reads(150, 150, 'I', Plain) + Reads(1, 2000, 'I', Plain) + Reads(150, 150, 'I', Plain)));
    ASSERT_TRUE(stray);
    EXPECT_EQ(stray->type, protal::ReadType::Single);
    // Nothing to read: no guess (a missing file, or a pipe, whose reads would be lost to the alignment).
    EXPECT_FALSE(Guess((dir.path / "missing.fq").string()));
    EXPECT_FALSE(Guess(dir.Write("empty.fq", "")));
}

TEST(ReadTypeDetection, LongReadsByTheirQuality) {
    ScratchDir dir;
    auto const ont = Guess(dir.Write("q18.fq", Reads(50, 4000, '3', Plain)));  // Q18
    ASSERT_TRUE(ont);
    EXPECT_EQ(ont->type, protal::ReadType::ONT);
    EXPECT_NE(ont->evidence.find("median read quality Q18.0"), std::string::npos) << ont->evidence;
    EXPECT_NE(ont->evidence.find("the first 50 reads up to 4.0 kb long, median 4.0 kb"), std::string::npos) << ont->evidence;
    auto const hifi = Guess(dir.Write("q35.fq", Reads(50, 4000, 'D', Plain)));  // Q35
    ASSERT_TRUE(hifi);
    EXPECT_EQ(hifi->type, protal::ReadType::PacBio);
    // A read's quality is its bases' mean error probability: a tenth of the bases at Q5 make Q40 reads Q15 (the mean of
    // the Phred values would be Q36.5).
    auto const mixed = Guess(dir.Write("mixed.fq", Reads(50, 4000, 'I', Plain, false, '&')));
    ASSERT_TRUE(mixed);
    EXPECT_EQ(mixed->type, protal::ReadType::ONT);
    EXPECT_NE(mixed->evidence.find("Q15.0"), std::string::npos) << mixed->evidence;
    // Without qualities: PacBio, and the evidence says why.
    auto const fasta = Guess(dir.Write("long.fa", Reads(50, 4000, 'I', Plain, true)));
    ASSERT_TRUE(fasta);
    EXPECT_EQ(fasta->type, protal::ReadType::PacBio);
    EXPECT_NE(fasta->evidence.find("without base qualities"), std::string::npos) << fasta->evidence;
    // Q0 at every base means unknown (pbsim3's PacBio reads): as without qualities.
    auto const q0 = Guess(dir.Write("q0.fq", Reads(50, 4000, '!', Plain)));
    ASSERT_TRUE(q0);
    EXPECT_EQ(q0->type, protal::ReadType::PacBio);
    EXPECT_NE(q0->evidence.find("(FASTA, or Q0 throughout)"), std::string::npos) << q0->evidence;
    // Gzipped reads read as plain ones.
    std::string const gz = (dir.path / "q18.fq.gz").string();
    {
        ogzstream os(gz.c_str());
        os << Reads(50, 4000, '3', Plain);
    }
    auto const zipped = Guess(gz);
    ASSERT_TRUE(zipped);
    EXPECT_EQ(zipped->type, protal::ReadType::ONT);
}

TEST(ReadTypeDetection, NamesBeforeQualities) {
    ScratchDir dir;
    // MinKNOW: a UUID, and runid=; dorado: a UUID. ONT reads, whatever their quality (duplex reads reach Q30).
    auto minknow = [](size_t i) {
        char id[64];
        std::snprintf(id, sizeof id, "%08zx-1c2d-4e5f-8a9b-0123456789ab runid=6f3a ch=12 start_time=2024-01-01T00:00:00Z", i);
        return std::string(id);
    };
    auto const ont = Guess(dir.Write("minknow.fq", Reads(40, 3000, 'D', minknow)));
    ASSERT_TRUE(ont);
    EXPECT_EQ(ont->type, protal::ReadType::ONT);
    EXPECT_NE(ont->evidence.find("named as MinKNOW and dorado name reads"), std::string::npos) << ont->evidence;
    // PacBio: movie/ZMW, also with one underscore in the movie name and more after the ZMW.
    auto const ccs = Guess(dir.Write("ccs.fq", Reads(40, 3000, '0', [](size_t i) { return "m64011_190830_220126/" + std::to_string(i) + "/ccs"; })));
    ASSERT_TRUE(ccs);
    EXPECT_EQ(ccs->type, protal::ReadType::PacBio);
    auto const revio = Guess(dir.Write("revio.fq", Reads(40, 3000, '0', [](size_t i) { return "m84011_220902_175841_s1/" + std::to_string(i) + "/ccs/fwd"; })));
    ASSERT_TRUE(revio);
    EXPECT_EQ(revio->type, protal::ReadType::PacBio);
    auto const one_underscore = Guess(dir.Write("m.fq", Reads(40, 3000, '0', [](size_t i) { return "m64001_000000/" + std::to_string(i) + "/ccs"; })));
    ASSERT_TRUE(one_underscore);
    EXPECT_EQ(one_underscore->type, protal::ReadType::PacBio);
    // Names like neither: by quality ('0' is Q15).
    auto const other = Guess(dir.Write("other.fq", Reads(40, 3000, '0', [](size_t i) { return "m_" + std::to_string(i); })));
    ASSERT_TRUE(other);
    EXPECT_EQ(other->type, protal::ReadType::ONT);
}
