// Unit tests for the features against false positives in complex communities (docs/claude/2026-10-06-false-positive-
// features): the unique k-mers of cores with one value, the read's crowding (ZN) and failed genes (ZF), the long reads'
// consensus tags (ZR), the database's species neighbours, and the features the profiler computes from them.
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include "Core/AlignmentStrategy.h"
#include "Hash/KmerLookup.h"
#include "Profiling/Profiler.h"
#include "SequenceUtils/SpeciesNeighbours.h"

using namespace protal;
namespace sn = protal::species_neighbours;

namespace {
    std::string RandomSequence(size_t length, std::mt19937& rng) {
        static constexpr char kBases[] = "ACGT";
        std::string seq(length, 'A');
        for (auto& c : seq) c = kBases[rng() % 4];
        return seq;
    }

    std::string Mutated(std::string seq, double rate, std::mt19937& rng) {
        static constexpr char kBases[] = "ACGT";
        std::bernoulli_distribution change(rate);
        for (auto& c : seq) {
            if (!change(rng)) continue;
            char other = c;
            while (other == c) other = kBases[rng() % 4];
            c = other;
        }
        return seq;
    }

    // A reference of taxa with genes 1..n each in a folder of its own, loaded.
    struct Reference {
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;
        std::map<uint32_t, std::vector<std::string>> genes;

        explicit Reference(std::map<uint32_t, std::vector<std::string>> taxa) : genes(std::move(taxa)) {
            dir = std::filesystem::temp_directory_path() / ("protal_cc_features_" + std::to_string(::getpid()) + "_" +
                                                            std::to_string(reinterpret_cast<uintptr_t>(this)));
            std::filesystem::create_directories(dir);
            std::ofstream fna(dir / "reference.fna");
            std::ofstream map(dir / "reference.map");
            size_t offset = 0;
            for (auto const& [taxid, seqs] : genes) {
                for (size_t i = 0; i < seqs.size(); i++) {
                    std::string const header = ">" + std::to_string(taxid) + "_" + std::to_string(i + 1) + "\n";
                    fna << header << seqs[i] << '\n';
                    map << taxid << '\t' << i + 1 << '\t' << offset + header.size() << '\t' << offset + header.size() + seqs[i].size() << '\n';
                    offset += header.size() + seqs[i].size() + 1;
                }
            }
            fna.close();
            map.close();
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
            loader->LoadAllGenomes();
        }

        std::string Header() const {
            std::string header = "@HD\tVN:1.6\n";
            for (auto const& [taxid, seqs] : genes) {
                for (size_t i = 0; i < seqs.size(); i++) {
                    header += "@SQ\tSN:" + std::to_string(taxid) + "_" + std::to_string(i + 1) + "\tLN:" + std::to_string(seqs[i].size()) + "\n";
                }
            }
            return header;
        }

        std::string Write(std::string const& name, std::string const& content) const {
            auto const path = (dir / name).string();
            std::ofstream(path) << content;
            return path;
        }

        ~Reference() { std::filesystem::remove_all(dir); }
    };

    constexpr int kPaired = 0x1, kBothAlign = 0x2, kReverse = 0x10, kRead1 = 0x40, kRead2 = 0x80, kSupplementary = 0x800;

    // A record of the gene's own bases from 0-based `start` (an exact match unless seq and cigar are given), with the
    // tags `tags` (tab-separated, without the leading tab; ZA:Z:* if empty).
    std::string Record(Reference const& ref, std::string const& qname, int flag, uint32_t taxid, uint32_t gene, int start, int length,
                       std::string const& tags = "", std::string seq = "", std::string cigar = "") {
        if (seq.empty()) seq = ref.genes.at(taxid)[gene - 1].substr(static_cast<size_t>(start), static_cast<size_t>(length));
        if (cigar.empty()) cigar = std::to_string(length) + "M";
        return qname + '\t' + std::to_string(flag) + '\t' + std::to_string(taxid) + "_" + std::to_string(gene) + '\t' + std::to_string(start + 1) +
               "\t60\t" + cigar + "\t*\t0\t0\t" + seq + '\t' + std::string(seq.size(), 'I') + '\t' + (tags.empty() ? "ZA:Z:*" : tags) + '\n';
    }

    std::map<std::string, double> Features(profiler::Taxon const& taxon) {
        std::map<std::string, double> features;
        for (auto const& [name, value] : profiler::TaxonFeatures(taxon)) features[name] = value;
        return features;
    }

    // Profiles a SAM (its text) as a run does, with taxa 1 and 2 in genus 10 and taxon 3 in genus 11.
    profiler::MicrobialProfile Profile(Reference const& ref, std::string const& sam, size_t threads = 1) {
        profiler::Profiler profiler(*ref.loader);
        profiler.SetDepthIdentityMargin(0.08);
        profiler::MicrobialProfile profile(*ref.loader);
        profile.SetGenera(std::make_shared<std::vector<uint32_t> const>(std::vector<uint32_t>{ 0, 10, 10, 11 }));
        auto const path = ref.Write("sample" + std::to_string(threads) + ".sam", sam);
        std::ostringstream rejected;
        std::string const error = profiler.ProfileSam(path, profile, std::optional<std::reference_wrapper<std::ostream>>{ rejected },
                                                      2, 2, 0.0, 15, 0, false, threads);
        EXPECT_EQ(error, "");
        return profile;
    }

    Reference ThreeTaxa(std::mt19937& rng) {
        std::map<uint32_t, std::vector<std::string>> taxa;
        for (uint32_t t : { 1u, 2u, 3u }) {
            for (int g = 0; g < 4; g++) taxa[t].push_back(RandomSequence(600, rng));
        }
        return Reference(std::move(taxa));
    }
}

TEST(ComplexCommunityFeatures, ASeedOfACoreWithOneValueKeepsItsFlagThroughTheSort) {
    LookupResult seed(7, 3, 120, 40, true, true, true);
    auto const back = LookupResult::FromSortKey(seed.SortKey());
    EXPECT_EQ(back.taxid, 7u);
    EXPECT_EQ(back.geneid, 3u);
    EXPECT_EQ(back.genepos, 120u);
    EXPECT_EQ(back.readpos, 40u);
    EXPECT_TRUE(back.unique);
    EXPECT_TRUE(back.unique_dist_two);
    EXPECT_TRUE(back.single);
    LookupResult plain(UINT32_MAX >> 1, (1u << 20) - 1, (1u << 30) + 5, 65000, true, false, false);
    auto const plain_back = LookupResult::FromSortKey(plain.SortKey());
    EXPECT_EQ(plain_back.taxid, UINT32_MAX >> 1);
    EXPECT_EQ(plain_back.geneid, (1u << 20) - 1);
    EXPECT_EQ(plain_back.genepos, (1u << 30) + 5);
    EXPECT_EQ(plain_back.readpos, 65000u);
    EXPECT_TRUE(plain_back.unique);
    EXPECT_FALSE(plain_back.unique_dist_two);
    EXPECT_FALSE(plain_back.single);
    // The flags sort last: taxon, gene, read and gene position first.
    EXPECT_LT(LookupResult(1, 1, 5, 2, true, true, true).SortKey(), LookupResult(1, 1, 6, 2).SortKey());
}

TEST(ComplexCommunityFeatures, FailedCandidatesKeepTheGeneOfTheirLongestAnchor) {
    using C = FailedCandidate;
    // Attempted in anchor order (longest first): taxon 3 on gene 7, then 1 on 2, 3 on 9, 2 on 4; taxon 2 aligned.
    EXPECT_EQ(FailedCandidates({ C(3, 7), C(1, 2), C(3, 9), C(2, 4) }, { 2 }), (std::vector<C>{ C(1, 2), C(3, 7) }));
    EXPECT_EQ(FailedCandidates({ 3, 1, 3, 2 }, { 2 }), (std::vector<C>{ 1, 3 }));  // genes unknown
    EXPECT_EQ(FailedTag({ C(1, 2), C(3, 7) }), "1:2,3:7");
    EXPECT_EQ(FailedTag({ 1, 3 }), "1,3");
    EXPECT_EQ(FailedTag({ C(1, 2), 3 }), "1:2,3");
    std::vector<std::pair<uint32_t, uint32_t>> seen;
    ForEachFailedCandidateGene("12:3,40,x,7:1,9:y,5:", [&](uint32_t t, uint32_t g) { seen.emplace_back(t, g); });
    EXPECT_EQ(seen, (std::vector<std::pair<uint32_t, uint32_t>>{ { 12, 3 }, { 40, 0 }, { 7, 1 }, { 5, 0 } }));
    std::vector<uint32_t> taxa;
    ForEachFailedCandidate("12:3,40", [&](uint32_t t) { taxa.push_back(t); });
    EXPECT_EQ(taxa, (std::vector<uint32_t>{ 12, 40 }));
}

TEST(ComplexCommunityFeatures, TheConsensusAndCrowdingTagsRoundTrip) {
    SamEntry sam;
    EXPECT_TRUE(UnmappedRecord(sam, "read1", { FailedCandidate(4, 2) }));
    EXPECT_EQ(sam.ToString(), "read1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tZU:i:0\tZT:i:0\tZF:Z:4:2");
    std::istringstream mapped("r\t0\t1_1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tZU:i:0\tZT:i:0\tZA:Z:*\tZF:Z:5:3\tZR:i:2\tZN:i:7\n"
                              "s\t0\t1_1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tZU:i:0\tZT:i:0\n");
    SamReader reader(mapped);
    SamEntry a, b;
    bool has_a = false, has_b = false;
    ASSERT_TRUE(reader.Next(a, b, has_a, has_b));
    ASSERT_TRUE(has_a);
    EXPECT_EQ(a.m_failed, "5:3");
    EXPECT_EQ(a.m_settled, 2);
    EXPECT_EQ(a.m_crowding, 7);
    EXPECT_EQ(a.ToString(), "r\t0\t1_1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tZU:i:0\tZT:i:0\tZA:Z:*\tZF:Z:5:3\tZR:i:2\tZN:i:7");
    ASSERT_TRUE(reader.Next(a, b, has_a, has_b));
    EXPECT_EQ(a.m_settled, 0);  // absent
    EXPECT_EQ(a.m_crowding, 0);
}

TEST(ComplexCommunityFeatures, CrowdingCountsTheTaxaWithAnAnchorNearTheLongest) {
    struct A {
        uint16_t total_length;
        uint32_t taxid;
    };
    std::vector<uint32_t> scratch;
    std::vector<A> anchors{ { 70, 3 }, { 100, 1 }, { 90, 2 }, { 85, 2 }, { 80, 4 } };  // 80 is 0.8 of 100
    EXPECT_EQ(SimpleAlignmentHandler::CrowdedTaxa(anchors, scratch), 3);
    EXPECT_EQ(SimpleAlignmentHandler::CrowdedTaxa(std::vector<A>{}, scratch), 0);
    EXPECT_EQ(SimpleAlignmentHandler::CrowdedTaxa(std::vector<A>{ { 50, 9 } }, scratch), 1);
}

TEST(ComplexCommunityFeatures, ACongenerThatFitsBetterThanItsDistanceAllowsIsUnexpected) {
    EXPECT_NEAR(profiler::PoissonCdf(0, 1.0), std::exp(-1.0), 1e-12);
    EXPECT_NEAR(profiler::PoissonCdf(2, 1.5), std::exp(-1.5) * (1 + 1.5 + 1.125), 1e-12);
    EXPECT_EQ(profiler::PoissonCdf(-1, 1.0), 0);
    EXPECT_EQ(profiler::PoissonCdf(3, 0.0), 1);
    EXPECT_EQ(profiler::PoissonCdf(5, 2000.0), 0);  // underflows to 0: the fit is as unexpected as can be
    // A read of 150 bases on a congener 5% away would differ ~7.5 bases more: a fit as good is unexpected.
    EXPECT_TRUE(profiler::UnexpectedFit(0, 0.05, 150));
    EXPECT_FALSE(profiler::UnexpectedFit(4, 0.05, 150));
    // On a congener 0.5% away a read fits with 0 or 1 more edit as a matter of course.
    EXPECT_FALSE(profiler::UnexpectedFit(0, 0.005, 150));
    EXPECT_FALSE(profiler::UnexpectedFit(1, 0.005, 150));
}

TEST(ComplexCommunityFeatures, TheSpeciesNeighboursTableKeepsTheNearestAndRoundTrips) {
    sn::Table table;
    EXPECT_TRUE(table.Empty());
    std::vector<sn::Neighbour> many;
    for (uint32_t i = 0; i < 20; i++) many.push_back({ 100 + i, 0.005f * static_cast<float>(i + 1) });  // 0.005 .. 0.1
    many.push_back({ 99, 0.2f });  // too far
    table.Set(1, many);
    table.Set(2, { { 1, 0.03f }, { 5, 0.012f } });
    table.Set(3, {});
    EXPECT_FALSE(table.Empty());
    EXPECT_EQ(table.Species(), 3u);
    ASSERT_EQ(table.Of(1).size(), sn::kMaxNeighbours);
    EXPECT_EQ(table.Of(1).front().taxid, 100u);
    EXPECT_EQ(table.Of(1).back().taxid, 115u);
    EXPECT_EQ(table.Of(2).front().taxid, 5u);  // nearest first
    EXPECT_TRUE(table.Of(3).empty());
    EXPECT_TRUE(table.Of(77).empty());
    EXPECT_TRUE(table.Has(3));
    EXPECT_FALSE(table.Has(4));
    EXPECT_EQ(table.Within(1, 0.0251), 5u);
    EXPECT_EQ(table.Within(2, 0.02), 1u);
    EXPECT_NEAR(table.Nearest(2), 0.012, 1e-6);
    EXPECT_EQ(table.Nearest(3), 1.0);
    // Listed: its distance; not listed: at least the farthest of a full list, kMaxDistance otherwise.
    EXPECT_NEAR(table.DistanceAtLeast(2, 1), 0.03, 1e-6);
    EXPECT_NEAR(table.DistanceAtLeast(1, 119), 0.08, 1e-6);
    EXPECT_NEAR(table.DistanceAtLeast(2, 9), sn::kMaxDistance, 1e-12);
    std::ostringstream os;
    table.Write(os);
    EXPECT_NE(os.str().find("taxid\tneighbours\n"), std::string::npos);
    EXPECT_NE(os.str().find("2\t5:0.012,1:0.03\n"), std::string::npos);
    EXPECT_NE(os.str().find("3\t\n"), std::string::npos);
    sn::Table back;
    std::istringstream is(os.str());
    ASSERT_EQ(back.Read(is), "");
    EXPECT_EQ(back.Species(), 3u);
    EXPECT_EQ(back.Pairs(), table.Pairs());
    for (uint32_t t : { 1u, 2u, 3u }) {
        ASSERT_EQ(back.Of(t).size(), table.Of(t).size());
        for (size_t i = 0; i < table.Of(t).size(); i++) EXPECT_EQ(back.Of(t)[i], table.Of(t)[i]);
    }
    sn::Table bad;
    std::istringstream bad_is("taxid\tneighbours\n2\t5-0.01\n");
    EXPECT_NE(bad.Read(bad_is), "");
}

TEST(ComplexCommunityFeatures, TheBuildsDistanceIsTheRunsDistance) {
    std::mt19937 rng(3);
    profiler::context::TaxonSketch a, b;
    for (uint32_t g = 1; g <= 14; g++) {
        auto const seq = RandomSequence(800, rng);
        a.emplace_back(g, profiler::context::GeneSketch(seq));
        if (g != 5) b.emplace_back(g, profiler::context::GeneSketch(Mutated(seq, 0.01 * (g % 4), rng)));
    }
    std::vector<double> scratch;
    EXPECT_EQ(profiler::context::SketchedTaxonDistance(a, b, scratch), profiler::context::SketchedTaxonDistances(a, b).distance);
    profiler::context::TaxonSketch few(a.begin(), a.begin() + 3);
    EXPECT_EQ(profiler::context::SketchedTaxonDistance(few, b, scratch), profiler::context::kFarDistance);
}

TEST(ComplexCommunityFeatures, GenesThatDivergeUnevenlyAreDispersed) {
    using G = profiler::GeneDivergence;
    // One divergence of 2% by the genes' factors (0.5 and 1.5) and no errors: the genes' differences as expected.
    std::vector<G> even{ { 10, 1000, 0, 0.5 }, { 30, 1000, 0, 1.5 }, { 20, 1000, 0, 1.0 }, { 10, 1000, 0, 0.5 } };
    EXPECT_NEAR(profiler::DivergenceDispersion(even), 0.0, 1e-9);
    // A mosaic: two genes identical, two diverged by 6%.
    std::vector<G> mosaic{ { 0, 1000, 0, 1 }, { 60, 1000, 0, 1 }, { 0, 1000, 0, 1 }, { 60, 1000, 0, 1 } };
    EXPECT_GT(profiler::DivergenceDispersion(mosaic), 10.0);
    // Fewer than kMinDispersionGenes genes with kMinGeneAligned bases: 1.
    std::vector<G> few{ { 0, 1000, 0, 1 }, { 60, 1000, 0, 1 }, { 5, 50, 0, 1 } };
    EXPECT_EQ(profiler::DivergenceDispersion(few), 1.0);
    // The errors are expected: reads of 1% errors and no divergence are not dispersed.
    std::vector<G> errors{ { 10, 1000, 10, 1 }, { 10, 1000, 10, 1 }, { 10, 1000, 10, 1 } };
    EXPECT_NEAR(profiler::DivergenceDispersion(errors), 0.0, 1e-9);
}

TEST(ComplexCommunityFeatures, GenesWithMoreFailedReadsThanRecordsAndComplementaryGenes) {
    tsl::robin_map<uint32_t, profiler::GeneRecords> records;
    records[1].records = 5;
    records[2].records = 3;
    tsl::robin_map<uint32_t, uint32_t> failed;
    EXPECT_EQ(profiler::FailedGeneShare(records, failed), 0);
    failed[2] = 4;  // more than its 3 records
    failed[3] = 2;  // a gene without records
    EXPECT_NEAR(profiler::FailedGeneShare(records, failed), 2.0 / 3, 1e-12);
    failed[2] = 3;
    EXPECT_NEAR(profiler::FailedGeneShare(records, failed), 1.0 / 3, 1e-12);

    std::vector<uint32_t> const genes{ 1, 2, 3, 4 };
    // Complementary: 0 shared hit genes where 1 is expected.
    EXPECT_NEAR(profiler::GeneOverlap({ 1, 2 }, { 3, 4 }, genes, genes), 0.5 / 1.5, 1e-12);
    // Both hit all genes: as expected.
    EXPECT_NEAR(profiler::GeneOverlap(genes, genes, genes, genes), 1.0, 1e-12);
    // Genes only one reference has do not count; no hit gene among the shared ones: 1.
    EXPECT_NEAR(profiler::GeneOverlap({ 1, 5 }, { 1 }, { 1, 2, 5 }, { 1, 2 }), (1 + 0.5) / (0.5 + 0.5), 1e-12);
    EXPECT_EQ(profiler::GeneOverlap({ 5 }, { 1 }, { 1, 2, 5 }, { 1, 2 }), 1.0);
}

TEST(ComplexCommunityFeatures, APairedSamplesFeaturesFromTheTagsAndTheMates) {
    std::mt19937 rng(11);
    Reference ref = ThreeTaxa(rng);
    sn::Table neighbours;
    neighbours.Set(1, { { 2, 0.048f } });
    neighbours.Set(2, { { 1, 0.048f } });
    neighbours.Set(3, {});
    ref.loader->SetSpeciesNeighbours(neighbours);
    std::string sam = ref.Header();
    // A: mates on taxa 1 and 2 (one genus): split for both. Mate 1 fits taxon 2 as well (0 more edits) where a read of
    // taxon 1 would differ by 0.048 x 100 = 4.8 bases (P(0) = 0.008): unexpected; it seeded on 4 taxa (ZN 4).
    sam += Record(ref, "A", kPaired | kRead1, 1, 1, 0, 100, "ZA:Z:2:0\tZN:i:4");
    sam += Record(ref, "A", kPaired | kRead2 | kReverse, 2, 1, 300, 100, "ZA:Z:*\tZN:i:1");
    // B: both mates on taxon 1's gene 2; mate 1 fits taxon 2 three edits worse: expected (P = 0.29).
    sam += Record(ref, "B", kPaired | kBothAlign | kRead1, 1, 2, 0, 100, "ZA:Z:2:3\tZN:i:2");
    sam += Record(ref, "B", kPaired | kBothAlign | kRead2 | kReverse, 1, 2, 300, 100, "ZA:Z:*");
    // C: mates on taxon 1 and taxon 3 (another genus): not split; the read failed on taxon 2's gene 3.
    sam += Record(ref, "C", kPaired | kRead1, 1, 3, 0, 100, "ZA:Z:*\tZF:Z:2:3");
    sam += Record(ref, "C", kPaired | kRead2 | kReverse, 3, 1, 300, 100);
    for (size_t threads : { 1, 3 }) {
        SCOPED_TRACE(threads);
        auto const profile = Profile(ref, sam, threads);
        auto const f1 = Features(profile.GetTaxa().at(1));
        auto const f2 = Features(profile.GetTaxa().at(2));
        auto const f3 = Features(profile.GetTaxa().at(3));
        EXPECT_NEAR(f1.at("split_fragment_share"), 1.0 / 3, 1e-12);
        EXPECT_EQ(f2.at("split_fragment_share"), 1.0);
        EXPECT_EQ(f3.at("split_fragment_share"), 0.0);
        // Taxon 1's four records: ZN 4, 2 and none twice: log2 4 = 2, log2 2 = 1.
        EXPECT_NEAR(f1.at("seed_crowding"), 1.5, 1e-9);
        EXPECT_EQ(f2.at("seed_crowding"), 0.0);  // log2 1
        EXPECT_NEAR(f1.at("unexpected_congener_fit_share"), 1.0 / 4, 1e-12);
        EXPECT_EQ(f2.at("unexpected_congener_fit_share"), 0.0);
        // Taxon 2: one record on gene 1, one failed read on gene 3.
        EXPECT_NEAR(f2.at("failed_gene_share"), 0.5, 1e-12);
        EXPECT_EQ(f1.at("failed_gene_share"), 0.0);
        // Taxon 1's partner is taxon 2 (its records' only congener alternative): of their 4 shared genes taxon 1 hits
        // 3 (1, 2, 3), taxon 2 one (1), both gene 1: as expected (3 x 1 / 4 = 0.75).
        EXPECT_NEAR(f1.at("congener_gene_overlap"), (1 + 0.5) / (0.75 + 0.5), 1e-12);
        EXPECT_EQ(f3.at("congener_gene_overlap"), 1.0);
        // The database around them: taxon 1's congener at 0.048.
        EXPECT_EQ(f1.at("db_congeners_01"), 0.0);
        EXPECT_EQ(f1.at("db_congeners_05"), 1.0);
        EXPECT_NEAR(f1.at("db_nearest_congener"), 0.048, 1e-6);
        EXPECT_EQ(f3.at("db_congeners_05"), 0.0);
        EXPECT_EQ(f3.at("db_nearest_congener"), 1.0);
        // Short reads have no ZR.
        EXPECT_EQ(f1.at("read_consensus_share"), 0.0);
        EXPECT_EQ(f1.at("read_inconsistent_share"), 0.0);
        // One read on a 600-base gene covers its 100 bases where 600 (1 - e^(-1/6)) are expected.
        EXPECT_NEAR(f3.at("breadth_ratio"), 100.0 / (600 * (1 - std::exp(-100.0 / 600))), 1e-9);
    }
    // Without the species neighbours: the database's features unknown, no fit unexpected.
    ref.loader->SetSpeciesNeighbours(sn::Table{});
    auto const f1 = Features(Profile(ref, sam).GetTaxa().at(1));
    EXPECT_EQ(f1.at("db_congeners_05"), sn::kUnknown);
    EXPECT_EQ(f1.at("db_nearest_congener"), sn::kUnknown);
    EXPECT_EQ(f1.at("unexpected_congener_fit_share"), 0.0);
}

TEST(ComplexCommunityFeatures, ALongReadsConsensusTags) {
    std::mt19937 rng(12);
    Reference ref = ThreeTaxa(rng);
    std::string sam = ref.Header();
    // One long read: genes 1 and 2 of taxon 1 (gene 1 settled by the read's consensus, ZR 1) and gene 3 of taxon 2,
    // clearly that taxon's (ZR 2).
    sam += Record(ref, "L", 0, 1, 1, 0, 400, "ZA:Z:*\tZR:i:1\tZN:i:2");
    sam += Record(ref, "L", kSupplementary, 1, 2, 0, 400);
    sam += Record(ref, "L", kSupplementary, 2, 3, 0, 400, "ZA:Z:*\tZR:i:2");
    auto const profile = Profile(ref, sam);
    auto const f1 = Features(profile.GetTaxa().at(1));
    auto const f2 = Features(profile.GetTaxa().at(2));
    EXPECT_NEAR(f1.at("read_consensus_share"), 0.5, 1e-12);
    EXPECT_EQ(f1.at("read_inconsistent_share"), 0.0);
    EXPECT_EQ(f2.at("read_inconsistent_share"), 1.0);
    EXPECT_EQ(f1.at("split_fragment_share"), 1.0);  // the read is split between the two congeners
    EXPECT_EQ(f2.at("split_fragment_share"), 1.0);
    EXPECT_NEAR(f1.at("seed_crowding"), 1.0, 1e-9);
}

TEST(ComplexCommunityFeatures, FixedDifferencesAndPolymorphicSitesOfTheConsensus) {
    std::mt19937 rng(13);
    Reference ref = ThreeTaxa(rng);
    auto const& gene = ref.genes.at(1)[0];
    auto mutate = [&](std::string seq, size_t pos) {
        seq[pos] = seq[pos] == 'A' ? 'C' : 'A';
        return seq;
    };
    std::string sam = ref.Header();
    // Four reads over bases 0-99 of taxon 1's gene 1, all with another base at 50 (a fixed difference) and two of them
    // at 70 too (a polymorphic site).
    for (int r = 0; r < 4; r++) {
        auto seq = mutate(gene.substr(0, 100), 50);
        std::string cigar = "50M1X49M";
        if (r < 2) {
            seq = mutate(seq, 70);
            cigar = "50M1X19M1X29M";
        }
        sam += Record(ref, "F" + std::to_string(r), 0, 1, 1, 0, 100, "", seq, cigar);
    }
    auto const f1 = Features(Profile(ref, sam).GetTaxa().at(1));
    EXPECT_NEAR(f1.at("fixed_difference_rate"), 1.0 / 100, 1e-12);
    EXPECT_NEAR(f1.at("polymorphic_site_rate"), 1.0 / 100, 1e-12);
}

TEST(ComplexCommunityFeatures, AnEmptyTaxonHasTheNeutralValues) {
    Genome no_genome(0);
    auto const f = Features(profiler::Taxon(no_genome));
    EXPECT_EQ(f.at("gene_divergence_dispersion"), 1.0);
    EXPECT_EQ(f.at("breadth_ratio"), 1.0);
    EXPECT_EQ(f.at("congener_gene_overlap"), 1.0);
    EXPECT_EQ(f.at("db_congeners_01"), sn::kUnknown);
    for (auto const* name : { "read_consensus_share", "read_inconsistent_share", "seed_crowding", "unexpected_congener_fit_share",
                              "split_fragment_share", "failed_gene_share", "fixed_difference_rate", "polymorphic_site_rate" }) {
        EXPECT_EQ(f.at(name), 0.0) << name;
    }
}
