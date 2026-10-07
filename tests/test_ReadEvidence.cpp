// Unit tests for the reads' evidence before the filters and what the database knows before any read
// (docs/claude/2026-10-03-false-positive-fixes): fragments before the MAPQ filter and the EM's fragments, the failed
// candidates (the ZF tag, on a read's first record or on an unmapped record, or the SAM header's line, read back and
// counted per taxon), and the species' priors (species_priors.tsv) as features.
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include "Profiling/Profiler.h"
#include "IO/AlignmentOutputHandler.h"
#include "SequenceUtils/SpeciesPriors.h"
#include "TestReference.h"

using namespace protal;
using protal::test::RandomSequence;
using Reference = protal::test::LoadedReference;

namespace {
    // A record of the gene's own bases from 0-based `start`; tags appended as given.
    std::string Record(Reference const& ref, std::string const& qname, int flag, uint32_t taxid, uint32_t gene, int start, int length,
                       int mapq = 60, std::string const& tags = "\tZA:Z:*") {
        std::string const seq = ref.genes.at(taxid)[gene - 1].substr(static_cast<size_t>(start), static_cast<size_t>(length));
        return qname + '\t' + std::to_string(flag) + '\t' + std::to_string(taxid) + "_" + std::to_string(gene) + '\t' + std::to_string(start + 1) +
               '\t' + std::to_string(mapq) + '\t' + std::to_string(length) + "M\t*\t0\t0\t" + seq + '\t' + std::string(seq.size(), 'I') + tags + '\n';
    }

    std::map<std::string, double> Features(profiler::Taxon const& taxon) {
        std::map<std::string, double> features;
        for (auto const& [name, value] : profiler::TaxonFeatures(taxon)) features[name] = value;
        return features;
    }

    // Profiles a SAM (its text) as a run does, on `threads` threads; on several, in chunks of about a read each.
    profiler::MicrobialProfile Profile(Reference const& ref, std::string const& sam, size_t threads = 1) {
        profiler::Profiler profiler(*ref.loader);
        profiler.SetDepthIdentityMargin(0.08);
        profiler.SetChunkBytes(200);
        profiler::MicrobialProfile profile(*ref.loader);
        auto const path = ref.Write("sample" + std::to_string(threads) + ".sam", sam);
        std::ostringstream rejected;
        std::string const error = profiler.ProfileSam(path, profile, std::optional<std::reference_wrapper<std::ostream>>{ rejected },
                                                      2, 2, 0.0, 15, 0, false, threads);
        EXPECT_EQ(error, "");
        return profile;
    }
}

TEST(ReadEvidence, FailedCandidatesAreTheAttemptedTaxaWithoutAnAlignment) {
    EXPECT_EQ(FailedCandidates({ 3, 1, 3, 2 }, { 2 }), (std::vector<FailedCandidate>{ 1, 3 }));
    EXPECT_EQ(FailedCandidates({ 2, 2 }, { 2, 5 }), (std::vector<FailedCandidate>{}));
    EXPECT_EQ(FailedCandidates({}, {}), (std::vector<FailedCandidate>{}));
    EXPECT_EQ(FailedTag({ 1, 3 }), "1,3");
    EXPECT_EQ(FailedTag({}), "");
    std::vector<uint32_t> seen;
    ForEachFailedCandidate("12,40,x,7", [&](uint32_t t) { seen.push_back(t); });
    EXPECT_EQ(seen, (std::vector<uint32_t>{ 12, 40, 7 }));
    seen.clear();
    ForEachFailedCandidate("*", [&](uint32_t t) { seen.push_back(t); });
    ForEachFailedCandidate("", [&](uint32_t t) { seen.push_back(t); });
    EXPECT_TRUE(seen.empty());
}

TEST(ReadEvidence, TheTagRoundTripsAndAnUnmappedRecordCarriesIt) {
    SamEntry sam;
    EXPECT_TRUE(UnmappedRecord(sam, "read7", { 4, 9 }));
    EXPECT_FALSE(UnmappedRecord(sam, "read8", {}));
    auto const line = sam.ToString();
    EXPECT_EQ(line, "read7\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:4,9");
    // Read back: the reader skips it as unmapped and counts its failed candidates per taxon.
    std::istringstream is("@HD\tVN:1.6\n" + line + "\nread9\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:9\n");
    SamReader reader(is);
    SamEntry a, b;
    bool has_a = false, has_b = false;
    EXPECT_FALSE(reader.Next(a, b, has_a, has_b));
    EXPECT_EQ(reader.Skipped().at("unmapped"), 2u);
    EXPECT_EQ(reader.FailedCandidates().at(4), 1u);
    EXPECT_EQ(reader.FailedCandidates().at(9), 2u);
    // A mapped record's tag is parsed into m_failed.
    std::istringstream mapped("r\t0\t1_1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tZU:i:0\tZT:i:0\tZA:Z:*\tZF:Z:5\n");
    SamReader reader2(mapped);
    ASSERT_TRUE(reader2.Next(a, b, has_a, has_b));
    ASSERT_TRUE(has_a);
    EXPECT_EQ(a.m_failed, "5");
    EXPECT_EQ(a.ToString().substr(a.ToString().size() - 6), "ZF:Z:5");
}

// The reader adds the SAM header's counts (FailedCandidatesLine, no line without failed candidates) and the unmapped
// records' tags up; a bad entry names its line. The output handler's line is tested in test_SamRoundTrip.cpp.
TEST(ReadEvidence, TheHeadersFailedCandidatesAddUpWithTheRecords) {
    EXPECT_EQ(FailedCandidatesLine({}), "");
    EXPECT_EQ(FailedCandidatesLine({ 0, 0 }), "");
    std::istringstream is("@HD\tVN:1.6\n" + FailedCandidatesLine({ 0, 0, 0, 0, 2, 0, 0, 0, 0, 1 }) +
                          "read9\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:9\n");
    SamReader reader(is);
    SamEntry a, b;
    bool has_a = false, has_b = false;
    EXPECT_FALSE(reader.Next(a, b, has_a, has_b));
    EXPECT_EQ(reader.FailedCandidates().at(4), 2u);
    EXPECT_EQ(reader.FailedCandidates().at(9), 2u);
    EXPECT_EQ(reader.Skipped().at("unmapped"), 1u);

    std::istringstream bad("@HD\tVN:1.6\n" + kSamFailedCandidatesComment + "4:2,9\n");
    SamReader bad_reader(bad);
    try {
        bad_reader.Next(a, b, has_a, has_b);
        ADD_FAILURE() << "a bad entry is not an error";
    } catch (SamFormatError const& e) {
        EXPECT_EQ(std::string(e.what()).rfind("line 2: ", 0), 0u) << e.what();
    }
}

// An unmapped record is skipped without a SamEntry (SamReader::SkipUnmapped), but one that cannot be parsed is the same
// error as before.
TEST(ReadEvidence, AnUnmappedRecordThatCannotBeParsedIsAnError) {
    for (std::string const line : { "r\t4\t*\tx\t0\t*\t*\t0\t0\t*\t*\n", "r\t4\t*\t0\t0\t*\t*\t0\n", "r\t4\t*\t0\t300\t*\t*\t0\t0\t*\t*\n" }) {
        std::istringstream is(line);
        SamReader reader(is);
        SamEntry a, b;
        bool has_a = false, has_b = false;
        EXPECT_THROW(reader.Next(a, b, has_a, has_b), SamFormatError) << line;
    }
}

// Two reads aligned nowhere: protal writes their failed candidates in the SAM header by default, as unmapped records
// with --write_unmapped_reads; either way they reach the taxa's failed_candidate_rate.
TEST(ReadEvidence, FragmentsBeforeTheFiltersTheEMsFragmentsAndTheFailedCandidates) {
    std::mt19937 rng(21);
    std::vector<std::string> one, two, three;
    for (int g = 0; g < 2; g++) {
        one.push_back(RandomSequence(500, rng));
        two.push_back(RandomSequence(500, rng));
        three.push_back(RandomSequence(500, rng));
    }
    Reference ref({ { 1, one }, { 2, two }, { 3, three } });
    std::string records;
    // Taxon 1: three reads kept, one more at MAPQ 0 (dropped by the filter but a best record), and one read that
    // names taxon 2 as a failed candidate on its first record.
    records += Record(ref, "a1", 0, 1, 1, 10, 100);
    records += Record(ref, "a2", 0, 1, 2, 10, 100);
    records += Record(ref, "a3", 0, 1, 1, 50, 100, 60, "\tZA:Z:*\tZF:Z:2");
    records += Record(ref, "a4", 0, 1, 2, 60, 100, 0);
    // Taxon 2: one kept read.
    records += Record(ref, "b1", 0, 2, 1, 10, 100);
    // Taxon 3: one read of its own and one that taxon 1 fits as well (no edit more). The EM gives that read the share
    // s = (1 + s) / ((1 + s) + (4 + 1 - s)) of it (taxon 3's weight against taxon 1's): s = 0.2, and taxon 3 keeps
    // 1.2 of its 2 reads.
    records += Record(ref, "c1", 0, 3, 1, 10, 100);
    records += Record(ref, "c2", 0, 3, 2, 10, 100, 60, "\tZA:Z:1:0");
    // Two reads that seeded on taxon 2 (one also on taxon 1) but aligned nowhere.
    std::string const unmapped = "u1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:2\n"
                                 "u2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:1,2\n";
    for (bool const in_header : { true, false }) {
        std::string const sam = in_header ? ref.Header() + FailedCandidatesLine({ 0, 1, 2 }) + records : ref.Header() + records + unmapped;
        for (size_t threads : { 1u, 3u }) {
            SCOPED_TRACE(std::string(in_header ? "in the header" : "unmapped records") + ", threads " + std::to_string(threads));
            auto const profile = Profile(ref, sam, threads);
            auto const f1 = Features(profile.GetTaxa().at(1));
            auto const f2 = Features(profile.GetTaxa().at(2));
            auto const f3 = Features(profile.GetTaxa().at(3));
            EXPECT_EQ(f1.at("fragments"), 3.0);
            EXPECT_EQ(f1.at("fragments_all"), 4.0);
            EXPECT_EQ(f1.at("em_own_share"), 1.0);  // no alternatives: the EM leaves it every read
            EXPECT_EQ(f1.at("em_fragments"), 4.0);
            EXPECT_NEAR(f3.at("em_own_share"), 0.6, 1e-6);
            EXPECT_NEAR(f3.at("em_kept_own_share"), 0.6, 1e-6);
            EXPECT_EQ(f3.at("fragments_all"), 2.0);
            EXPECT_NEAR(f3.at("em_fragments"), 1.2, 2e-6);
            // Taxon 1 failed for u2 only: 1 of 1 + 4.
            EXPECT_NEAR(f1.at("failed_candidate_rate"), 1.0 / 5, 1e-12);
            // Taxon 2 failed for a3, u1 and u2: 3 of 3 + 1.
            EXPECT_EQ(f2.at("fragments_all"), 1.0);
            EXPECT_NEAR(f2.at("failed_candidate_rate"), 3.0 / 4, 1e-12);
            EXPECT_EQ(f3.at("failed_candidate_rate"), 0.0);
        }
    }
}

TEST(ReadEvidence, SpeciesPriorsAreReadAndBecomeFeatures) {
    species_priors::Table table;
    std::istringstream is("taxid\trep_genome\tmarkers\tduplicate_markers\tcheckm_completeness\tcheckm_contamination\tani_radius\tmean_intra_ani\tmin_intra_ani\tclustered_genomes\n"
                          "1\tGCF_1\t120\t3\t99.5\t1.2\t95\t98.7\t96.1\t250\n"
                          "2\tGCA_2\t53\t0\t-1\tnone\t\t-1\t-1\t1\n");
    ASSERT_EQ(table.Read(is), "");
    EXPECT_EQ(table.Size(), 2u);
    EXPECT_EQ(table.Informative(), 2u);
    EXPECT_EQ(table.Get(1).duplicate_markers, 3.0);
    EXPECT_EQ(table.Get(2).completeness, species_priors::kUnknown);
    EXPECT_EQ(table.Get(2).clustered_genomes, 1.0);
    EXPECT_EQ(table.Get(7).markers, species_priors::kUnknown);
    std::istringstream bad("1\t2\t3\n");
    EXPECT_NE(species_priors::Table().Read(bad), "");
    std::istringstream bad_number("1\tG\tx\t0\t0\t0\t0\t0\t0\t0\n");
    EXPECT_NE(species_priors::Table().Read(bad_number), "");

    std::mt19937 rng(22);
    Reference ref({ { 1, { RandomSequence(400, rng) } }, { 2, { RandomSequence(400, rng) } } });
    ref.loader->SetSpeciesPriors(table);
    std::string sam = ref.Header();
    sam += Record(ref, "a", 0, 1, 1, 10, 100);
    sam += Record(ref, "b", 0, 2, 1, 10, 100);
    auto const profile = Profile(ref, sam);
    auto const f1 = Features(profile.GetTaxa().at(1));
    auto const f2 = Features(profile.GetTaxa().at(2));
    EXPECT_NEAR(f1.at("rep_duplicate_share"), 3.0 / 120, 1e-12);
    EXPECT_EQ(f1.at("rep_completeness"), 99.5);
    EXPECT_EQ(f1.at("rep_contamination"), 1.2);
    EXPECT_EQ(f1.at("cluster_ani_radius"), 95.0);
    EXPECT_EQ(f1.at("cluster_min_ani"), 96.1);
    EXPECT_NEAR(f1.at("cluster_genomes_log10"), std::log10(250.0), 1e-12);
    EXPECT_EQ(f2.at("rep_duplicate_share"), 0.0);
    EXPECT_EQ(f2.at("rep_completeness"), species_priors::kUnknown);
    EXPECT_EQ(f2.at("cluster_ani_radius"), species_priors::kUnknown);
    EXPECT_EQ(f2.at("cluster_genomes_log10"), 0.0);
    // Without the table every value is unknown.
    Reference plain({ { 1, { RandomSequence(400, rng) } } });
    auto const none = Profile(plain, plain.Header() + Record(plain, "a", 0, 1, 1, 10, 100));
    EXPECT_EQ(Features(none.GetTaxa().at(1)).at("rep_duplicate_share"), species_priors::kUnknown);
}
