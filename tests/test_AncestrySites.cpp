// Unit tests of AncestrySites.h: the sites where a species' gene copy differs from its congener's, found along their
// shared 12-mers; what a record covers of them from its CIGAR; and the per-run cache over the species neighbours.
#include <gtest/gtest.h>
#include <random>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>
#include "SequenceUtils/AncestrySites.h"

using namespace protal;
namespace an = protal::ancestry;
namespace sn = protal::species_neighbours;

namespace {
    std::string RandomSequence(size_t length, std::mt19937& rng) {
        static constexpr char kBases[] = "ACGT";
        std::string seq(length, 'A');
        for (auto& c : seq) c = kBases[rng() % 4];
        return seq;
    }

    char Other(char c, std::mt19937& rng) {
        static constexpr char kBases[] = "ACGT";
        char other = c;
        while (other == c) other = kBases[rng() % 4];
        return other;
    }

    // own with substitutions at the given positions; returns the new bases.
    std::string Substituted(std::string seq, std::vector<size_t> const& positions, std::mt19937& rng) {
        for (auto const p : positions) seq[p] = Other(seq[p], rng);
        return seq;
    }
}

TEST(AncestrySites, CompareFindsTheSubstitutions) {
    std::mt19937 rng(11);
    auto const own = RandomSequence(900, rng);
    std::vector<size_t> changed;
    for (size_t p = 30; p < 870; p += 29) changed.push_back(p);
    auto const other = Substituted(own, changed, rng);
    auto const sites = an::Compare(own, other);
    ASSERT_TRUE(sites.has_value());
    EXPECT_EQ(std::vector<size_t>(sites->positions.begin(), sites->positions.end()), changed);
    for (size_t i = 0; i < changed.size(); i++) EXPECT_EQ(sites->bases[i], an::Code(other[changed[i]]));
    EXPECT_GT(sites->compared, 800u);
    EXPECT_NEAR(sites->identity, 1.0 - static_cast<double>(changed.size()) / static_cast<double>(sites->compared), 1e-6);
    EXPECT_EQ(sites->Count(0, 900), changed.size());
    EXPECT_EQ(sites->Count(30, 31), 1u);
    EXPECT_EQ(sites->Count(31, 59), 0u);
    EXPECT_EQ(sites->Find(59), 1u);
    EXPECT_EQ(sites->Find(60), SIZE_MAX);
    // Identical copies: no site, everything compared.
    auto const same = an::Compare(own, own);
    ASSERT_TRUE(same.has_value());
    EXPECT_TRUE(same->Empty());
    EXPECT_EQ(same->compared, 900u);
}

TEST(AncestrySites, CompareStopsAtAnIndelAndRejectsUnrelatedCopies) {
    std::mt19937 rng(12);
    auto const own = RandomSequence(900, rng);
    std::vector<size_t> changed = { 100, 200, 300, 400 };
    auto other = Substituted(own, changed, rng);
    other.insert(500, "ACGTACG");  // an insertion: the stretch after it leaves the main diagonal
    auto const sites = an::Compare(own, other);
    ASSERT_TRUE(sites.has_value());
    EXPECT_EQ(std::vector<size_t>(sites->positions.begin(), sites->positions.end()), changed);
    EXPECT_LT(sites->compared, 600u);  // the first half only
    EXPECT_GT(sites->compared, 450u);
    // Two unrelated sequences share no unique 12-mer.
    EXPECT_FALSE(an::Compare(own, RandomSequence(900, rng)).has_value());
    // Too short to pair, or too long to index.
    EXPECT_FALSE(an::Compare("ACGTACGTAC", "ACGTACGTAC").has_value());
    EXPECT_FALSE(an::Compare(std::string(70000, 'A'), std::string(70000, 'A')).has_value());
}

TEST(AncestrySites, CountWalksTheCigar) {
    an::Sites sites;
    sites.positions = { 10, 20, 30, 45 };
    sites.bases = { an::Code('T'), an::Code('G'), an::Code('C'), an::Code('A') };
    sites.congener = 7;
    // A record from position 1 (1-based) covering 50 bases with no mismatch: every site covered, all the species' base.
    auto c = an::Count(sites, "50M", 1, std::string(50, 'A'));
    EXPECT_EQ(c.sites, 4u);
    EXPECT_EQ(c.agree, 4u);
    EXPECT_EQ(c.congener, 0u);
    // A mismatch at site 20 with the congener's base (G), one at 30 with another base (T: neither), one off-site.
    std::string seq(50, 'A');
    seq[20] = 'G';
    seq[30] = 'T';
    seq[35] = 'C';
    c = an::Count(sites, "20M1X9M1X4M1X14M", 1, seq);
    EXPECT_EQ(c.sites, 4u);
    EXPECT_EQ(c.agree, 2u);
    EXPECT_EQ(c.congener, 1u);
    // From position 11: sites 10 and 20 covered by 11M (reference 10-20), 30 under the deletion (21-30, skipped), 45
    // in the last run (31-60).
    c = an::Count(sites, "11M10D30M", 11, std::string(41, 'A'));
    EXPECT_EQ(c.sites, 3u);
    EXPECT_EQ(c.agree, 3u);
    // Soft clips and insertions consume the read only: from position 8 (reference 7), 3S then 3M (7-9), the X at
    // reference 10 is the read's base 6 (the congener's T), the 43M (11-53) cover the other three sites.
    seq = std::string(50, 'A');
    seq[6] = 'T';
    c = an::Count(sites, "3S3M1X43M", 8, seq);
    EXPECT_EQ(c.sites, 4u);
    EXPECT_EQ(c.agree, 3u);
    EXPECT_EQ(c.congener, 1u);
    // An insertion before the X: the read's base shifts, the reference position does not.
    seq = std::string(52, 'A');
    seq[8] = 'T';
    c = an::Count(sites, "3S3M2I1X43M", 8, seq);
    EXPECT_EQ(c.sites, 4u);
    EXPECT_EQ(c.congener, 1u);
    // No sequence, or no sites: nothing.
    EXPECT_EQ(an::Count(sites, "50M", 1, "*").sites, 0u);
    EXPECT_EQ(an::Count(an::Sites{}, "50M", 1, std::string(50, 'A')).sites, 0u);
}

namespace {
    // A stand-in for the GenomeLoader: the copies of a few species.
    struct FakeSequence {
        std::string seq;
        std::string_view View() const { return seq; }
    };
    struct FakeGene {
        std::string seq;
        FakeSequence Sequence() const { return { seq }; }
    };
    struct FakeGenomes {
        std::unordered_map<uint64_t, FakeGene> genes;
        void Set(uint32_t taxid, uint32_t gene, std::string seq) { genes[(static_cast<uint64_t>(taxid) << 32) | gene] = { std::move(seq) }; }
        bool HasGene(uint32_t taxid, uint32_t gene) const { return genes.contains((static_cast<uint64_t>(taxid) << 32) | gene); }
        FakeGene const& GetGeneOMP(uint32_t taxid, uint32_t gene) { return genes.at((static_cast<uint64_t>(taxid) << 32) | gene); }
    };
}

TEST(AncestrySites, CacheTakesTheNearestCongenerWithThePairingCopy) {
    std::mt19937 rng(13);
    auto const own = RandomSequence(800, rng);
    FakeGenomes genomes;
    genomes.Set(1, 5, own);
    genomes.Set(3, 5, Substituted(own, { 50, 150, 250 }, rng));  // the second-nearest congener: the only one with the gene
    genomes.Set(4, 5, Substituted(own, { 60, 160, 260, 360 }, rng));
    genomes.Set(1, 6, own);
    genomes.Set(2, 6, RandomSequence(800, rng));  // the nearest congener's copy of gene 6 does not pair
    genomes.Set(3, 6, Substituted(own, { 70 }, rng));
    sn::Table neighbours;
    neighbours.Set(1, { { 2, 0.02f }, { 3, 0.04f }, { 4, 0.06f } });
    an::Cache cache;
    auto const s5 = cache.Get(1, 5, genomes, neighbours);
    EXPECT_EQ(s5->congener, 3u);
    EXPECT_EQ(s5->positions.size(), 3u);
    EXPECT_EQ(cache.Get(1, 5, genomes, neighbours).get(), s5.get());  // computed once
    auto const s6 = cache.Get(1, 6, genomes, neighbours);
    EXPECT_EQ(s6->congener, 3u);
    EXPECT_EQ(s6->positions.size(), 1u);
    // A species without neighbours, or without the gene: empty; without the table: empty and nothing cached.
    EXPECT_TRUE(cache.Get(3, 5, genomes, neighbours)->Empty());
    EXPECT_TRUE(cache.Get(1, 7, genomes, neighbours)->Empty());
    EXPECT_EQ(cache.Size(), 4u);
    EXPECT_TRUE(cache.Get(1, 5, genomes, sn::Table{})->Empty());
    cache.Clear();
    EXPECT_EQ(cache.Size(), 0u);
}
