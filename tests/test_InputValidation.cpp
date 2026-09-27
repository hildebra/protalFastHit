// Unit tests for the database and input files protal reads: malformed files must stop with a clear
// message (exit 8) instead of being half-read, and the index must know which reference it was built for.
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <unistd.h>
#include "Hash/Seedmap.h"
#include "Profiling/Profiler.h"
#include "SequenceUtils/GenomeLoader.h"
#include "SequenceUtils/KmerUtils.h"
#include "Taxonomy/Taxonomy.h"
#include "Utilities/ReferenceFingerprint.h"

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
