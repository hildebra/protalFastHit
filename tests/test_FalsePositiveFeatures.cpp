// Unit tests for what was added against false positives (docs/claude/2026-10-03-false-positive-anatomy): the sample's
// depth as a feature, the reads' divergence by gene conservation and by codon position, lost mates, and the suspect
// gene copies (GeneIncongruence.h) that --build finds and a run leaves out.
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <unistd.h>
#include "Profiling/Profiler.h"
#include "SequenceUtils/GeneIncongruence.h"
#include "TestReference.h"

using namespace protal;
namespace gi = protal::gene_incongruence;
using protal::test::Mutated;
using protal::test::RandomSequence;
using Reference = protal::test::LoadedReference;

namespace {
    constexpr int kPaired = 0x1, kBothAlign = 0x2, kMateUnmapped = 0x8, kReverse = 0x10, kMateReverse = 0x20, kRead1 = 0x40, kRead2 = 0x80;

    // A record of the gene's own bases from 0-based `start` (an exact match unless a CIGAR is given).
    std::string Record(Reference const& ref, std::string const& qname, int flag, uint32_t taxid, uint32_t gene, int start, int length,
                       int mapq = 60, std::string cigar = "", std::string seq = "") {
        if (seq.empty()) seq = ref.genes.at(taxid)[gene - 1].substr(static_cast<size_t>(start), static_cast<size_t>(length));
        if (cigar.empty()) cigar = std::to_string(length) + "M";
        return qname + '\t' + std::to_string(flag) + '\t' + std::to_string(taxid) + "_" + std::to_string(gene) + '\t' + std::to_string(start + 1) +
               '\t' + std::to_string(mapq) + '\t' + cigar + "\t*\t0\t0\t" + seq + '\t' + std::string(seq.size(), 'I') + "\tZA:Z:*\n";
    }

    SamEntry MakeSam(uint32_t taxid, uint32_t geneid, std::string seq, POS_t pos, std::string cigar, int flag = 0) {
        SamEntry sam;
        sam.m_qname = "read";
        sam.m_flag = static_cast<FLAG_t>(flag);
        sam.m_rname = std::to_string(taxid) + "_" + std::to_string(geneid);
        sam.m_pos = pos;
        sam.m_mapq = 60;
        sam.m_cigar = std::move(cigar);
        sam.m_qual = std::string(seq.size(), 'I');
        sam.m_seq = std::move(seq);
        sam.m_alternatives = "*";
        return sam;
    }

    std::map<std::string, double> Features(profiler::Taxon const& taxon) {
        std::map<std::string, double> features;
        for (auto const& [name, value] : profiler::TaxonFeatures(taxon)) features[name] = value;
        return features;
    }

    // Profiles a SAM (its text) with the profiler, as a run does, on `threads` threads; on several, in chunks of about
    // a read each.
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

TEST(FalsePositiveFeatures, MismatchesAreCountedByCodonPositionOfTheReference) {
    // Position 1 (1-based) is codon position 1; X at 0-based reference positions 2, 5, 8 are third positions.
    EXPECT_EQ(profiler::MismatchesByCodonPosition("9M", 1), (std::pair<uint64_t, uint64_t>{ 0, 0 }));
    EXPECT_EQ(profiler::MismatchesByCodonPosition("2M1X2M1X2M1X", 1), (std::pair<uint64_t, uint64_t>{ 3, 3 }));
    EXPECT_EQ(profiler::MismatchesByCodonPosition("1X1X1X", 1), (std::pair<uint64_t, uint64_t>{ 3, 1 }));
    // From position 3 (0-based 2, a third position): X at 2 (third), then after 1M at 4 (second).
    EXPECT_EQ(profiler::MismatchesByCodonPosition("1X1M1X", 3), (std::pair<uint64_t, uint64_t>{ 2, 1 }));
    // Insertions and clips consume no reference, deletions do: X after 1M2D is at 0-based 3 (first position),
    // after 1M2I at 1 (second).
    EXPECT_EQ(profiler::MismatchesByCodonPosition("1M2D1X", 1), (std::pair<uint64_t, uint64_t>{ 1, 0 }));
    EXPECT_EQ(profiler::MismatchesByCodonPosition("2S1M2I1X", 1), (std::pair<uint64_t, uint64_t>{ 1, 0 }));
    EXPECT_EQ(profiler::MismatchesByCodonPosition("1M1X", 2), (std::pair<uint64_t, uint64_t>{ 1, 1 }));
}

TEST(FalsePositiveFeatures, TheRoomForAMateAndAFragmentsSpan) {
    SamEntry forward = MakeSam(1, 1, std::string(100, 'A'), 11, "100M");
    EXPECT_EQ(profiler::MateRoom(forward, 1000), 990u);  // from 0-based 10 to the gene's end
    EXPECT_EQ(profiler::MateRoom(forward, 0), 0u);
    SamEntry reverse = MakeSam(1, 1, std::string(100, 'A'), 801, "100M", kReverse);
    EXPECT_EQ(profiler::MateRoom(reverse, 1000), 900u);  // from the gene's start to the record's end
    SamEntry with_deletion = MakeSam(1, 1, std::string(100, 'A'), 801, "50M5D50M", kReverse);
    EXPECT_EQ(profiler::MateRoom(with_deletion, 1000), 905u);
    SamEntry mate = MakeSam(1, 1, std::string(100, 'A'), 251, "100M", kReverse);
    EXPECT_EQ(profiler::FragmentSpan(forward, mate), 340u);
    EXPECT_EQ(profiler::FragmentSpan(mate, forward), 340u);
}

TEST(FalsePositiveFeatures, TheMeanErrorProbabilityComesFromTheQualities) {
    SamEntry sam = MakeSam(1, 1, "ACGT", 1, "4M");
    sam.m_qual = "IIII";  // Q40
    EXPECT_NEAR(*profiler::MeanErrorProbability(sam), 1e-4, 1e-12);
    sam.m_qual = "+!+!";  // Q10 and Q0
    EXPECT_NEAR(*profiler::MeanErrorProbability(sam), (0.1 + 1 + 0.1 + 1) / 4, 1e-12);
    sam.m_qual = "*";
    EXPECT_FALSE(profiler::MeanErrorProbability(sam).has_value());
    // ReadExcess is the divergence beyond it: the differences per aligned base less the mean error probability of the
    // record's bases; none without qualities or aligned bases.
    EXPECT_FALSE(profiler::ReadExcess(sam).has_value());
    sam.m_qual = "IIII";
    sam.m_cigar = "3M1X";
    EXPECT_NEAR(*profiler::ReadExcess(sam), 0.25 - 1e-4, 1e-6);
    sam.m_cigar = "10M";
    sam.m_qual = std::string(10, 'I');  // Q40: 0.0001 per base
    EXPECT_NEAR(*profiler::ReadExcess(sam), -1e-4, 1e-7);
    sam.m_cigar = "8M2X";
    sam.m_qual = std::string(10, '+');  // Q10: 0.1
    EXPECT_NEAR(*profiler::ReadExcess(sam), 0.1, 1e-6);
    sam.m_cigar = "2S4M1D4M";  // a difference in 9 aligned bases; the clipped bases' qualities count too
    sam.m_qual = std::string(10, '5');  // Q20: 0.01
    EXPECT_NEAR(*profiler::ReadExcess(sam), 1.0 / 9 - 0.01, 1e-6);
    sam.m_cigar = "10S";
    EXPECT_FALSE(profiler::ReadExcess(sam).has_value());
}

TEST(FalsePositiveFeatures, AProfilesTaxaCarryTheSamplesDepthAndTheirReadsDivergenceByConservationAndCodonPosition) {
    std::mt19937 rng(7);
    std::vector<std::string> one, two;
    for (int g = 0; g < 4; g++) {
        one.push_back(RandomSequence(600, rng));
        two.push_back(RandomSequence(600, rng));
    }
    Reference ref({ { 1, one }, { 2, two } });
    // Genes 1 and 2 conserved (factor 0.5), 3 and 4 fast (1.5).
    gene_conservation::Table factors;
    for (uint64_t g = 1; g <= 4; g++) factors.Set(g, g <= 2 ? 0.5 : 1.5, 3);
    ref.loader->SetGeneConservation(factors);

    profiler::MicrobialProfile profile(*ref.loader);
    int read = 0;
    // Taxon 1: 9 reads, 100 bases each, on every gene: 2 mismatches per read on the conserved genes (at third codon
    // positions), 6 on the fast ones (at first positions) at Q40: excesses 0.02 and 0.06, scaled 0.04 each.
    for (int r = 0; r < 9; r++) {
        uint32_t const gene = static_cast<uint32_t>(r % 4) + 1;
        bool const conserved = gene <= 2;
        std::string const cigar = conserved ? "2M1X47M1X49M" : "1X2M1X2M1X2M1X2M1X2M1X84M";
        auto const exact = MakeSam(1, gene, one[gene - 1].substr(0, 100), 1, "100M");
        auto noted = exact;  // what the records' evidence sees; the profile takes the read as the gene's own bases
        noted.m_cigar = cigar;
        profile.NoteRecord(1, gene, noted);
        ASSERT_TRUE(profile.AddSam(1, static_cast<int>(gene), exact, 0.98, true, read++));
    }
    // Taxon 2: one exact read.
    auto sam = MakeSam(2, 1, two[0].substr(0, 100), 1, "100M");
    profile.NoteRecord(2, 1, sam);
    ASSERT_TRUE(profile.AddSam(2, 1, sam, 1.0, true, read++));
    profile.ApplyRecordEvidence();

    auto const f1 = Features(profile.GetTaxa().at(1));
    auto const f2 = Features(profile.GetTaxa().at(2));
    // The sample's depth: log10 of 10 fragments, the same for every taxon.
    EXPECT_NEAR(f1.at("sample_log_fragments"), 1.0, 1e-12);
    EXPECT_NEAR(f2.at("sample_log_fragments"), 1.0, 1e-12);
    EXPECT_EQ(profile.Fragments(), 10u);
    // The sample's complexity, the same for both: two taxa, the low-identity share of their fragments, and the median
    // identity over both (neither has 10 fragments), which is taxon 1's: it holds 9 of the 10 fragments.
    auto const& t1 = profile.GetTaxa().at(1);
    auto const& t2 = profile.GetTaxa().at(2);
    for (auto const* f : { &f1, &f2 }) {
        EXPECT_NEAR(f->at("sample_log_taxa"), std::log10(2.0), 1e-12);
        EXPECT_NEAR(f->at("sample_low_identity"), (9 * t1.LowIdentityShare() + t2.LowIdentityShare()) / 10, 1e-12);
        EXPECT_EQ(f->at("sample_identity"), t1.BaseIdentity());
    }
    // Taxon 1's excess: medians over 9 reads (5 conserved at 0.02, 4 fast at 0.06) -> 0.02 - 1e-4; scaled all 0.04 - e.
    EXPECT_NEAR(f1.at("excess_median"), 0.02 - 1e-4, 1e-6);
    EXPECT_NEAR(f1.at("excess_scaled_median"), (0.02 - 1e-4) / 0.5, 1e-6);
    // Conserved genes diverge at 0.02 beyond errors, fast ones at 0.06: log2 of their ratio, each + 0.001 (a default
    // feature: the number is pinned), log2(0.0209 / 0.0609).
    EXPECT_NEAR(f1.at("excess_conserved_fast_ratio"), -1.54294, 1e-4);
    // Mismatches: 5 reads x 2 at third positions, 4 reads x 6 at first positions.
    EXPECT_NEAR(f1.at("third_position_share"), 10.0 / 34, 1e-12);
    // Taxon 2 has no mismatch: the neutral third, and nothing beyond its errors.
    EXPECT_NEAR(f2.at("third_position_share"), 1.0 / 3, 1e-12);
    EXPECT_NEAR(f2.at("excess_scaled_median"), -1e-4 / 0.5, 1e-9);
    EXPECT_EQ(f2.at("excess_conserved_fast_ratio"), 0.0);  // too few bases on either kind
    // No paired reads went through the SAM path here: no mates judged.
    EXPECT_EQ(f1.at("mate_lost_share"), 0.0);
}

TEST(FalsePositiveFeatures, LostMatesAreCountedWhereTheFragmentWouldHaveFitTheGene) {
    std::mt19937 rng(8);
    std::vector<std::string> one, two;
    for (int g = 0; g < 2; g++) {
        one.push_back(RandomSequence(1200, rng));
        two.push_back(RandomSequence(1200, rng));
    }
    Reference ref({ { 1, one }, { 2, two } });
    std::string sam = ref.Header();
    // Taxon 1: 60 pairs on gene 1 with both mates (spans 300-360: the sample's fragment length), and 10 single mates,
    // forward: 5 at the start of the gene (room about 1180 for the mate: lost), 5 at its end (room about 160: not
    // judged).
    for (int p = 0; p < 60; p++) {
        int const start = 10 + p;
        sam += Record(ref, "pair" + std::to_string(p), kPaired | kBothAlign | kRead1 | kMateReverse, 1, 1, start, 100);
        sam += Record(ref, "pair" + std::to_string(p), kPaired | kBothAlign | kRead2 | kReverse, 1, 1, start + 200 + p, 100);
    }
    for (int s = 0; s < 5; s++) sam += Record(ref, "lost" + std::to_string(s), kPaired | kRead1 | kMateUnmapped, 1, 1, 20 + s, 100);
    for (int s = 0; s < 5; s++) sam += Record(ref, "edge" + std::to_string(s), kPaired | kRead1 | kMateUnmapped, 1, 1, 1040 + s, 100);
    // Taxon 2: a pair whose mates land on two taxa (read 2 on taxon 1's gene 2): taxon 2's mate has room, taxon 1's too.
    sam += Record(ref, "split", kPaired | kRead1 | kMateReverse, 2, 1, 10, 100);
    sam += Record(ref, "split", kPaired | kRead2 | kReverse, 1, 2, 900, 100);
    // A single-end read on taxon 2: not judged.
    sam += Record(ref, "single", 0, 2, 2, 10, 100);

    for (size_t threads : { 1u, 3u }) {
        SCOPED_TRACE(threads);
        auto const profile = Profile(ref, sam, threads);
        auto const f1 = Features(profile.GetTaxa().at(1));
        auto const f2 = Features(profile.GetTaxa().at(2));
        // Taxon 1: 60 linked fragments, 5 lost (room about 1180 >= the fragments' 95th percentile ~357), 5 unjudged, and the
        // split pair's mate on gene 2 (room 1000, reverse: lost too): 6 of 66.
        EXPECT_NEAR(f1.at("mate_lost_share"), 6.0 / 66, 1e-12);
        // Taxon 2: the split pair's read 1 (room 1190, lost); the single-end read is not judged.
        EXPECT_NEAR(f2.at("mate_lost_share"), 1.0, 1e-12);
        EXPECT_EQ(profile.GetTaxa().at(1).Fragments(), 71u);
        EXPECT_EQ(profile.GetTaxa().at(2).Fragments(), 2u);
    }
}

TEST(FalsePositiveFeatures, RecordsOnSuspectCopiesAreLeftOutAsIfTheReadHadNotAligned) {
    std::mt19937 rng(9);
    std::vector<std::string> one, two;
    for (int g = 0; g < 3; g++) {
        one.push_back(RandomSequence(500, rng));
        two.push_back(RandomSequence(500, rng));
    }
    Reference ref({ { 1, one }, { 2, two } });
    std::string sam = ref.Header();
    for (int r = 0; r < 6; r++) sam += Record(ref, "a" + std::to_string(r), 0, 1, static_cast<uint32_t>(r % 3) + 1, 10 + r, 100);
    sam += Record(ref, "b", 0, 2, 2, 10, 100);  // taxon 2's only read, on its suspect copy of gene 2
    sam += Record(ref, "c", kPaired | kBothAlign | kRead1 | kMateReverse, 1, 3, 10, 100);  // a pair with one mate on a suspect copy
    sam += Record(ref, "c", kPaired | kBothAlign | kRead2 | kReverse, 1, 3, 200, 100);

    auto const plain = Profile(ref, sam);
    EXPECT_EQ(plain.GetTaxa().count(2), 1u);
    EXPECT_EQ(plain.GetTaxa().at(1).Fragments(), 7u);
    EXPECT_EQ(plain.SuspectRecords(), 0u);

    gi::Table table;
    table.Add({ 2, 2, 1, gi::Rank::Family, 0.005f, 2.0f });
    table.Add({ 1, 3, 2, gi::Rank::Order, 0.0f, 2.0f });
    ref.loader->SetSuspectCopies(table);
    EXPECT_TRUE(ref.loader->IsSuspectCopy(2, 2));
    EXPECT_FALSE(ref.loader->IsSuspectCopy(2, 1));
    for (size_t threads : { 1u, 2u }) {
        SCOPED_TRACE(threads);
        auto const profile = Profile(ref, sam, threads);
        // Taxon 2 has no record left, taxon 1 lost its reads on gene 3 (two of the six, and the pair): 4 fragments.
        EXPECT_EQ(profile.GetTaxa().count(2), 0u);
        EXPECT_EQ(profile.GetTaxa().at(1).Fragments(), 4u);
        EXPECT_EQ(profile.GetTaxa().at(1).GetGenes().count(3), 0u);
        EXPECT_EQ(profile.SuspectRecords(), 5u);  // b, two a's on gene 3, c's two mates
        EXPECT_NEAR(Features(profile.GetTaxa().at(1)).at("sample_log_fragments"), std::log10(4.0), 1e-12);
    }
}

TEST(GeneIncongruence, LineagesComeFromTheTaxonomyFile) {
    std::string const dmp = "id\tparent_id\texternal_id\tname\trank\tlevel\trep_genome\n"
                            "1\t10\t0\ts__A one\tspecies\t7\tG1\n"
                            "2\t10\t0\ts__A two\tspecies\t7\tG2\n"
                            "3\t11\t0\ts__B one\tspecies\t7\tG3\n"
                            "10\t20\t0\tg__A\tgenus\t6\t\n"
                            "11\t21\t0\tg__B\tgenus\t6\t\n"
                            "20\t30\t0\tf__X\tfamily\t5\t\n"
                            "21\t31\t0\tf__Y\tfamily\t5\t\n"
                            "30\t40\t0\to__O\torder\t4\t\n"
                            "31\t41\t0\to__P\torder\t4\t\n"
                            "40\t50\t0\tc__C\tclass\t3\t\n"
                            "41\t50\t0\tc__D\tclass\t3\t\n"
                            "50\t60\t0\tp__P\tphylum\t2\t\n"
                            "60\t60\t0\td__Bacteria\tdomain\t1\t\n";
    std::istringstream is(dmp);
    std::unordered_map<uint32_t, gi::Lineage> lineages;
    ASSERT_EQ(gi::LineagesFromTaxonomy(is, lineages), "");
    EXPECT_EQ(lineages.at(1).At(gi::Rank::Genus), 10u);
    EXPECT_EQ(lineages.at(1).At(gi::Rank::Family), 20u);
    EXPECT_EQ(lineages.at(1).At(gi::Rank::Phylum), 50u);
    EXPECT_EQ(lineages.at(1).At(gi::Rank::Domain), 60u);
    EXPECT_EQ(gi::SharedRank(lineages.at(1), lineages.at(2)), gi::Rank::Genus);
    EXPECT_EQ(gi::SharedRank(lineages.at(1), lineages.at(3)), gi::Rank::Phylum);
    EXPECT_EQ(gi::SharedRank(lineages.at(1), gi::Lineage{}), gi::Rank::None);
    EXPECT_EQ(gi::RankName(gi::Rank::Family), "family");
    EXPECT_EQ(*gi::RankFromName("order"), gi::Rank::Order);
    EXPECT_FALSE(gi::RankFromName("kingdom").has_value());
    std::istringstream empty("");
    EXPECT_NE(gi::LineagesFromTaxonomy(empty, lineages), "");
}

TEST(GeneIncongruence, ACopyNearAnotherGenusAndFarFromItsCongenersIsSuspect) {
    std::mt19937 rng(10);
    // Genus 10: species 1, 2, 3; genus 11 (same family): species 4, 5; genus 12 (other phylum): species 6.
    std::unordered_map<uint32_t, gi::Lineage> lineages;
    auto lineage = [](uint32_t genus, uint32_t family, uint32_t order, uint32_t cls, uint32_t phylum) {
        gi::Lineage l;
        l.ancestor = { genus, family, order, cls, phylum, 60 };
        return l;
    };
    for (uint32_t t : { 1u, 2u, 3u }) lineages[t] = lineage(10, 20, 30, 40, 50);
    for (uint32_t t : { 4u, 5u }) lineages[t] = lineage(11, 20, 30, 40, 50);
    lineages[6] = lineage(12, 22, 32, 42, 52);
    auto sketch = [](std::string const& seq) { return gene_conservation::BottomSketch(seq, gi::kSketchSize); };
    std::vector<std::vector<gi::Copy>> copies(4);
    // Gene 1: each genus has its own gene, species within a genus 2% apart: nothing suspect.
    auto const g1a = RandomSequence(900, rng), g1b = RandomSequence(900, rng), g1c = RandomSequence(900, rng);
    copies[1] = { { 1, sketch(g1a) }, { 2, sketch(Mutated(g1a, 0.02, rng)) }, { 3, sketch(Mutated(g1a, 0.02, rng)) },
                  { 4, sketch(g1b) }, { 5, sketch(Mutated(g1b, 0.02, rng)) }, { 6, sketch(g1c) } };
    // Gene 2: species 3's copy is species 4's (another genus of the family) while its congeners' are unrelated:
    // suspect; species 4's congener 5 is at 1%, so 4 is not.
    auto const g2a = RandomSequence(900, rng), g2b = RandomSequence(900, rng);
    copies[2] = { { 1, sketch(g2a) }, { 2, sketch(Mutated(g2a, 0.02, rng)) }, { 3, sketch(Mutated(g2b, 0.003, rng)) },
                  { 4, sketch(g2b) }, { 5, sketch(Mutated(g2b, 0.01, rng)) }, { 6, sketch(Mutated(g2a, 0.08, rng)) } };
    // Gene 3: species 6 (alone in its genus, another phylum) shares its copy with species 1, whose congeners are at
    // 3.5%: 6 lies inside genus 10's cluster, tighter than the cluster, so it is suspect; and so is 1, nearer to 6 than
    // to its own congeners (which of the two carries the other's copy, the gene cannot tell). Species 4 (its genus has
    // no other copy of gene 3 either) is unrelated.
    auto const g3a = RandomSequence(900, rng);
    copies[3] = { { 1, sketch(g3a) }, { 2, sketch(Mutated(g3a, 0.035, rng)) }, { 3, sketch(Mutated(g3a, 0.035, rng)) },
                  { 4, sketch(RandomSequence(900, rng)) }, { 6, sketch(Mutated(g3a, 0.002, rng)) } };
    // Gene 4: a slow gene a young family shares: species 6, alone in its genus, at 1.5% from species 1, whose
    // congeners are at 1%: not inside the cluster by the margin, so not suspect; nor is 1.
    auto const g4a = RandomSequence(900, rng);
    copies.push_back({ { 1, sketch(g4a) }, { 2, sketch(Mutated(g4a, 0.01, rng)) }, { 3, sketch(Mutated(g4a, 0.01, rng)) },
                       { 6, sketch(Mutated(g4a, 0.015, rng)) } });
    for (int threads : { 1, 3 }) {
        SCOPED_TRACE(threads);
        auto const result = gi::Scan(copies, lineages, gi::kDefaultSuspectDistance, threads);
        EXPECT_EQ(result.copies, 21u);
        EXPECT_EQ(result.genes, 4u);
        ASSERT_EQ(result.suspects.size(), 3u);
        EXPECT_EQ(result.suspects[0].taxid, 3u);
        EXPECT_EQ(result.suspects[0].geneid, 2u);
        EXPECT_EQ(result.suspects[0].partner, 4u);
        EXPECT_EQ(result.suspects[0].rank, gi::Rank::Family);
        EXPECT_LT(result.suspects[0].distance, 0.01);
        EXPECT_GT(result.suspects[0].congener_distance, 0.05);
        EXPECT_EQ(result.suspects[1].taxid, 1u);
        EXPECT_EQ(result.suspects[1].geneid, 3u);
        EXPECT_EQ(result.suspects[1].partner, 6u);
        EXPECT_NEAR(result.suspects[1].congener_distance, 0.035, 0.01);
        EXPECT_EQ(result.suspects[2].taxid, 6u);
        EXPECT_EQ(result.suspects[2].geneid, 3u);
        EXPECT_EQ(result.suspects[2].partner, 1u);
        EXPECT_EQ(result.suspects[2].rank, gi::Rank::Domain);
        EXPECT_GE(result.suspects[2].congener_distance, 2.0f);  // no congener
        // Pairs across genera within kReportDistance: gene 2 (3,4), (3,5); gene 3 (1,6), (2,6), (3,6); gene 4 (1,6),
        // (2,6), (3,6).
        EXPECT_EQ(result.pairs.size(), 8u);
        EXPECT_TRUE(std::is_sorted(result.pairs.begin(), result.pairs.end(), [](gi::Pair const& a, gi::Pair const& b) {
            return std::tie(a.geneid, a.taxid_a, a.taxid_b) < std::tie(b.geneid, b.taxid_a, b.taxid_b);
        }));
        // A distance no copy reaches finds none.
        EXPECT_EQ(gi::Scan(copies, lineages, -1.0, threads).suspects.size(), 0u);
    }
}

namespace {
    // gene_incongruence::Scan before it skipped candidates that cannot be within kReportDistance (c9b46f9): every
    // candidate merged in full, genes in parallel. The new scan must give exactly what this gives.
    gi::Result ScanEveryCandidate(std::vector<std::vector<gi::Copy>> const& copies_by_gene,
                                  std::unordered_map<uint32_t, gi::Lineage> const& lineages, double suspect_distance) {
        gi::Result result;
        gi::Lineage const no_lineage;
        auto lineage_of = [&](uint32_t taxid) -> gi::Lineage const& {
            auto const it = lineages.find(taxid);
            return it == lineages.end() ? no_lineage : it->second;
        };
        for (size_t geneid = 0; geneid < copies_by_gene.size(); geneid++) {
            auto const& copies_here = copies_by_gene[geneid];
            size_t const n = copies_here.size();
            result.copies += n;
            if (n < 2) continue;
            result.genes++;
            std::vector<std::pair<uint32_t, uint32_t>> entries;
            for (size_t i = 0; i < n; i++) {
                for (uint32_t const h : copies_here[i].sketch) entries.emplace_back(h, static_cast<uint32_t>(i));
            }
            std::sort(entries.begin(), entries.end());
            std::unordered_map<uint32_t, size_t> per_genus;
            std::vector<uint32_t> genus_of(n, 0);
            for (size_t i = 0; i < n; i++) {
                genus_of[i] = lineage_of(copies_here[i].taxid).At(gi::Rank::Genus);
                if (genus_of[i] != 0) per_genus[genus_of[i]]++;
            }
            std::vector<float> congener(n, 2), foreign(n, 2);
            std::vector<size_t> partner(n, n);
            std::vector<gi::Rank> foreign_rank(n, gi::Rank::None);
            std::vector<uint32_t> found;
            for (size_t i = 0; i < n; i++) {
                auto const& copy = copies_here[i];
                found.clear();
                size_t const probes = std::min(gi::kProbeHashes, copy.sketch.size());
                for (size_t p = 0; p < probes; p++) {
                    auto const [lo, hi] = std::equal_range(entries.begin(), entries.end(), std::make_pair(copy.sketch[p], uint32_t{0}),
                        [](auto const& a, auto const& b) { return a.first < b.first; });
                    if (static_cast<size_t>(hi - lo) > gi::kMaxBucket) continue;
                    for (auto it = lo; it != hi; ++it) {
                        if (it->second != i) found.push_back(it->second);
                    }
                }
                std::sort(found.begin(), found.end());
                found.erase(std::unique(found.begin(), found.end()), found.end());
                result.candidates += found.size();
                auto const& lineage = lineage_of(copy.taxid);
                for (uint32_t const j : found) {
                    auto const& other = copies_here[j];
                    if (other.taxid == copy.taxid) continue;
                    auto const d = static_cast<float>(gene_conservation::SketchDistance(copy.sketch, other.sketch));
                    if (d > gi::kReportDistance) continue;
                    gi::Rank const rank = gi::SharedRank(lineage, lineage_of(other.taxid));
                    if (rank == gi::Rank::Genus) {
                        congener[i] = std::min(congener[i], d);
                        continue;
                    }
                    if (d < foreign[i] || (d == foreign[i] && other.taxid < copies_here[partner[i]].taxid)) {
                        foreign[i] = d;
                        partner[i] = j;
                        foreign_rank[i] = rank;
                    }
                    if (copy.taxid < other.taxid) result.pairs.push_back({ static_cast<uint32_t>(geneid), copy.taxid, other.taxid, rank, d });
                    else result.pairs.push_back({ static_cast<uint32_t>(geneid), other.taxid, copy.taxid, rank, d });
                }
            }
            for (size_t i = 0; i < n; i++) {
                if (foreign[i] > suspect_distance || partner[i] >= n) continue;
                bool const has_congeners = genus_of[i] != 0 && per_genus[genus_of[i]] > 1;
                bool suspect;
                if (has_congeners) {
                    suspect = congener[i] > foreign[i] + gi::kCongenerMargin;
                } else {
                    size_t const p = partner[i];
                    bool const partner_has_congeners = genus_of[p] != 0 && per_genus[genus_of[p]] > 1;
                    suspect = partner_has_congeners && congener[p] < 2 && foreign[i] + gi::kCongenerMargin <= congener[p];
                }
                if (suspect) {
                    result.suspects.push_back({ copies_here[i].taxid, static_cast<uint32_t>(geneid), copies_here[partner[i]].taxid,
                                                foreign_rank[i], foreign[i], congener[i] });
                }
            }
        }
        gi::Finish(result);
        return result;
    }
}

TEST(GeneIncongruence, SkippingCandidatesThatCannotBeNearChangesNothing) {
    // 40 genera of 1-12 species in 8 families; per gene, each genus a variant of the gene's family sequence, its species
    // 0-6% apart, a few copies taken from another genus (transfers), short genes (sketches below 128 hashes) and genes of
    // repeats: suspects and near pairs at every distance, and candidates far beyond kReportDistance.
    std::mt19937 rng(77);
    std::unordered_map<uint32_t, gi::Lineage> lineages;
    std::vector<std::pair<uint32_t, uint32_t>> species;  // taxid, genus
    uint32_t taxid = 1;
    for (uint32_t genus = 0; genus < 40; genus++) {
        uint32_t const family = genus % 8;
        size_t const members = 1 + rng() % 12;
        for (size_t m = 0; m < members; m++, taxid++) {
            gi::Lineage l;
            l.ancestor = { 1000 + genus, 2000 + family, 3000 + family / 2, 4000 + family / 4, 5000, 6000 };
            lineages[taxid] = l;
            species.emplace_back(taxid, genus);
        }
    }
    std::uniform_real_distribution<double> unit(0, 1);
    std::vector<std::vector<gi::Copy>> copies(7);
    for (size_t gene = 1; gene < copies.size(); gene++) {
        size_t const length = gene == 3 ? 90 : gene == 5 ? 400 : 1200;  // gene 3: fewer than 128 distinct 12-mers
        std::vector<std::string> family_seq(8), genus_seq(40);
        for (auto& f : family_seq) f = RandomSequence(length, rng);
        if (gene == 6) {  // repeats: a short motif many times
            std::string const motif = RandomSequence(30, rng);
            for (auto& f : family_seq) { f.clear(); while (f.size() < length) f += motif; f = Mutated(f, 0.02, rng); }
        }
        for (uint32_t g = 0; g < 40; g++) genus_seq[g] = Mutated(family_seq[g % 8], 0.03 + 0.1 * unit(rng), rng);
        for (auto const& [t, g] : species) {
            if (unit(rng) < 0.05) continue;  // the gene missing in this species
            uint32_t const source = unit(rng) < 0.04 ? static_cast<uint32_t>(rng() % 40) : g;  // a transfer, now and then
            copies[gene].push_back({ t, gene_conservation::BottomSketch(Mutated(genus_seq[source], 0.06 * unit(rng), rng), gi::kSketchSize) });
        }
        std::shuffle(copies[gene].begin(), copies[gene].end(), rng);
    }
    auto const expected = ScanEveryCandidate(copies, lineages, gi::kDefaultSuspectDistance);
    ASSERT_GT(expected.suspects.size(), 0u);
    ASSERT_GT(expected.pairs.size(), 0u);
    for (int threads : { 1, 4 }) {
        SCOPED_TRACE(threads);
        auto const result = gi::Scan(copies, lineages, gi::kDefaultSuspectDistance, threads);
        EXPECT_EQ(result.copies, expected.copies);
        EXPECT_EQ(result.genes, expected.genes);
        EXPECT_EQ(result.candidates, expected.candidates);
        EXPECT_LT(result.compared, result.candidates);  // some skipped
        ASSERT_EQ(result.suspects.size(), expected.suspects.size());
        for (size_t i = 0; i < result.suspects.size(); i++) {
            auto const& a = result.suspects[i];
            auto const& b = expected.suspects[i];
            EXPECT_EQ(std::tie(a.taxid, a.geneid, a.partner, a.rank, a.distance, a.congener_distance),
                      std::tie(b.taxid, b.geneid, b.partner, b.rank, b.distance, b.congener_distance)) << i;
        }
        ASSERT_EQ(result.pairs.size(), expected.pairs.size());
        for (size_t i = 0; i < result.pairs.size(); i++) {
            auto const& a = result.pairs[i];
            auto const& b = expected.pairs[i];
            EXPECT_EQ(std::tie(a.geneid, a.taxid_a, a.taxid_b, a.rank, a.distance), std::tie(b.geneid, b.taxid_a, b.taxid_b, b.rank, b.distance)) << i;
        }
    }
}

TEST(GeneIncongruence, TheBoundOnSharedHashesHoldsForEveryPair) {
    // SharedBound never falls below the hashes SketchDistance's merge finds shared, also when several of a sketch's
    // hashes fall on one bit of its signature; MergeShared gives SketchDistance's count, or nothing below `need`.
    std::mt19937 rng(5);
    auto signature = [](std::vector<uint32_t> const& sketch, std::vector<uint64_t>& sig) {
        sig.assign(gi::detail::kSignatureWords, 0);
        for (uint32_t const h : sketch) sig[(h & 1023u) >> 6] |= uint64_t{1} << (h & 63u);
        uint32_t bits = 0;
        for (uint64_t const w : sig) bits += static_cast<uint32_t>(__builtin_popcountll(w));
        return static_cast<uint32_t>(sketch.size()) - bits;
    };
    auto const need = gi::detail::SharedNeeded(gi::kSketchSize);
    EXPECT_EQ(need[0], 1u);
    for (int trial = 0; trial < 2000; trial++) {
        std::string const base = RandomSequence(300 + rng() % 900, rng);
        auto const a = gene_conservation::BottomSketch(base, gi::kSketchSize);
        auto const b = gene_conservation::BottomSketch(Mutated(base, 0.002 * (rng() % 60), rng), gi::kSketchSize);
        std::vector<uint64_t> sa, sb;
        uint32_t const ea = signature(a, sa), eb = signature(b, sb);
        size_t const s = std::min(a.size(), b.size());
        size_t const full = *gi::detail::MergeShared(a.data(), a.size(), b.data(), b.size(), s, 0);
        EXPECT_GE(gi::detail::SharedBound(sa.data(), sb.data(), ea, eb), full);
        EXPECT_EQ(gene_conservation::SketchDistanceOf(full, s), gene_conservation::SketchDistance(a, b));
        auto const stopped = gi::detail::MergeShared(a.data(), a.size(), b.data(), b.size(), s, need[s]);
        if (full >= need[s]) EXPECT_EQ(stopped, std::optional<size_t>(full));
        else EXPECT_FALSE(stopped.has_value());
        // A distance within kReportDistance is never below need.
        if (!(static_cast<float>(gene_conservation::SketchDistance(a, b)) > gi::kReportDistance)) EXPECT_GE(full, need[s]);
    }
}

TEST(GeneIncongruence, TheTableRoundTripsAndReportsItsSpecies) {
    gi::Table table;
    table.Add({ 3, 2, 4, gi::Rank::Family, 0.0049f, 2.0f });
    table.Add({ 6, 3, 1, gi::Rank::Domain, 0.0f, 0.021f });
    table.Add({ 6, 3, 1, gi::Rank::Domain, 0.0f, 0.021f });  // once
    EXPECT_EQ(table.Size(), 2u);
    EXPECT_EQ(table.Species(), 2u);
    EXPECT_TRUE(table.Contains(3, 2));
    EXPECT_FALSE(table.Contains(3, 3));
    std::ostringstream os;
    table.Write(os);
    EXPECT_NE(os.str().find("taxid\tgeneid\tpartner_taxid\tshared_rank\tdistance\tcongener_distance\n"), std::string::npos);
    EXPECT_NE(os.str().find("3\t2\t4\tfamily\t0.0049\tnone\n"), std::string::npos);
    gi::Table back;
    std::istringstream is(os.str());
    ASSERT_EQ(back.Read(is), "");
    EXPECT_EQ(back.Size(), 2u);
    EXPECT_TRUE(back.Contains(6, 3));
    EXPECT_EQ(back.Rows()[1].rank, gi::Rank::Domain);
    EXPECT_NEAR(back.Rows()[1].congener_distance, 0.021, 1e-6);
    std::istringstream bad("1\t2\t3\n");
    EXPECT_NE(gi::Table().Read(bad), "");
    std::istringstream bad_rank("1\t2\t3\tkingdom\t0\t0\n");
    EXPECT_NE(gi::Table().Read(bad_rank), "");
    // The report lists pairs with their verdicts.
    gi::Result result;
    result.pairs.push_back({ 2, 3, 4, gi::Rank::Family, 0.0049f });
    result.suspects.push_back({ 3, 2, 4, gi::Rank::Family, 0.0049f, 2.0f });
    std::ostringstream report;
    gi::WriteReport(report, result);
    EXPECT_NE(report.str().find("2\t3\t4\tfamily\t0.0049\t1\t0\n"), std::string::npos);
}
