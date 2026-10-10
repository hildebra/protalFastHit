// Unit tests of AncestrySites.h: the sites where a species' gene copy differs from its congener's, found along their
// shared 12-mers; the consensus over several congeners; what a record covers of them from its CIGAR; and the per-run
// cache over the species neighbours and the congener gaps.
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

    // The k-th (0, 1, 2) base other than c.
    char Alt(char c, size_t k) {
        std::string others;
        for (char const b : std::string_view("ACGT")) {
            if (b != c) others += b;
        }
        return others[k];
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
    EXPECT_EQ(sites->congeners, 1u);
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
    // The comparison says which positions were compared.
    auto const comparison = an::CompareCopies(own, other);
    ASSERT_TRUE(comparison.has_value());
    EXPECT_EQ(comparison->covered.size(), 900u);
    EXPECT_EQ(static_cast<size_t>(std::count(comparison->covered.begin(), comparison->covered.end(), 1)), comparison->sites.compared);
}

// The chain follows an indel onto the new diagonal: the stretch past it is compared, and the indel is a site of the
// comparison when it is three bases or more.
TEST(AncestrySites, CompareChainsPastAnIndelAndRecordsIt) {
    std::mt19937 rng(12);
    auto const own = RandomSequence(900, rng);
    std::vector<size_t> changed = { 100, 200, 300, 400, 700 };
    auto other = Substituted(own, changed, rng);
    other.insert(500, "ACGTACG");  // 7 bases the congener has and the species lacks
    auto sites = an::Compare(own, other);
    ASSERT_TRUE(sites.has_value());
    EXPECT_EQ(std::vector<size_t>(sites->positions.begin(), sites->positions.end()), changed);  // 700 lies past it
    EXPECT_GT(sites->compared, 850u);
    ASSERT_EQ(sites->indels.size(), 1u);
    EXPECT_EQ(sites->indels[0].length, -7);
    EXPECT_GE(sites->indels[0].position, 495u);
    EXPECT_LE(sites->indels[0].position, 505u);
    // The species' extra bases: 9 the congener lacks at 320.
    auto shorter = Substituted(own, { 150, 650 }, rng);
    shorter.erase(320, 9);
    sites = an::Compare(own, shorter);
    ASSERT_TRUE(sites.has_value());
    EXPECT_EQ(sites->positions, (std::vector<uint16_t>{ 150, 650 }));
    EXPECT_GT(sites->compared, 850u);
    ASSERT_EQ(sites->indels.size(), 1u);
    EXPECT_EQ(sites->indels[0].length, 9);
    EXPECT_GE(sites->indels[0].position, 315u);
    EXPECT_LE(sites->indels[0].position, 325u);
    // A single base lost: chained, but no indel site.
    auto frameshift = Substituted(own, { 100, 800 }, rng);
    frameshift.erase(450, 1);
    sites = an::Compare(own, frameshift);
    ASSERT_TRUE(sites.has_value());
    EXPECT_EQ(sites->positions, (std::vector<uint16_t>{ 100, 800 }));
    EXPECT_GT(sites->compared, 850u);
    EXPECT_TRUE(sites->indels.empty());
    // An insertion beyond kMaxIndel ends the chain: the first half only, as before.
    auto far = own;
    far.insert(500, RandomSequence(70, rng));
    sites = an::Compare(own, far);
    ASSERT_TRUE(sites.has_value());
    EXPECT_LT(sites->compared, 600u);
    EXPECT_GT(sites->compared, 450u);
    EXPECT_TRUE(sites->indels.empty());
    // Two unrelated sequences share no unique 12-mer.
    EXPECT_FALSE(an::Compare(own, RandomSequence(900, rng)).has_value());
    // Too short to pair, or too long to index.
    EXPECT_FALSE(an::Compare("ACGTACGTAC", "ACGTACGTAC").has_value());
    EXPECT_FALSE(an::Compare(std::string(70000, 'A'), std::string(70000, 'A')).has_value());
}

// The consensus: at a position compared with three congeners or more, a site where nine in ten of them carry one base
// other than the species' (every one of three or four); with fewer compared, the nearest congener's difference.
TEST(AncestrySites, ConsensusKeepsTheSitesTheCongenersShare) {
    std::mt19937 rng(14);
    auto const own = RandomSequence(900, rng);
    auto set = [&own](std::string& seq, size_t p, size_t k) { seq[p] = Alt(own[p], k); };
    // A: the nearest (three differences, all of it compared); B and C end at an insertion of 70 bases at 860, D at 700
    // (beyond kMaxIndel: the chain ends there).
    std::string a = own, b = own, c = own, d = own;
    for (auto* s : { &a, &b, &c, &d }) set(*s, 100, 0);          // all four: a site
    for (auto* s : { &b, &c, &d }) set(*s, 200, 0);              // three of four: none (0.75)
    set(b, 600, 0);                                              // three of four, two bases: none
    set(c, 600, 1);
    set(d, 600, 1);
    for (auto* s : { &a, &b, &c }) set(*s, 800, 2);              // the three compared there (D ended): a site
    set(b, 850, 0);                                              // one of three: none
    set(a, 880, 1);                                              // A alone compared there: its difference, a site
    set(b, 890, 1);                                              // B is not compared there
    b.insert(860, RandomSequence(70, rng));
    c.insert(860, RandomSequence(70, rng));
    d.insert(700, RandomSequence(70, rng));
    std::vector<an::Comparison> comparisons;
    for (auto const& [taxid, seq] : { std::pair{ 11u, &a }, std::pair{ 12u, &b }, std::pair{ 13u, &c }, std::pair{ 14u, &d } }) {
        auto comparison = an::CompareCopies(own, *seq);
        ASSERT_TRUE(comparison.has_value()) << taxid;
        comparison->sites.congener = taxid;
        comparisons.push_back(std::move(*comparison));
    }
    auto const sites = an::Consensus(own.size(), comparisons);
    EXPECT_EQ(sites.positions, (std::vector<uint16_t>{ 100, 800, 880 }));
    EXPECT_EQ(sites.bases, (std::vector<uint8_t>{ an::Code(Alt(own[100], 0)), an::Code(Alt(own[800], 2)), an::Code(Alt(own[880], 1)) }));
    EXPECT_EQ(sites.congener, 11u);
    EXPECT_EQ(sites.congeners, 4u);
    EXPECT_GT(sites.compared, 850u);
    EXPECT_NEAR(sites.identity, comparisons[0].sites.identity, 1e-7);
    // With B alone: its differences where compared, as 0.7.9 had (850 lies between that substitution and B's insertion
    // at 860, so no paired 12-mer reaches it); without any: nothing.
    auto const alone = an::Consensus(own.size(), { comparisons[1] });
    EXPECT_EQ(alone.positions, (std::vector<uint16_t>{ 100, 200, 600, 800 }));
    EXPECT_EQ(alone.congeners, 1u);
    EXPECT_TRUE(an::Consensus(own.size(), {}).Empty());
}

// From six congeners compared, one may carry a third base at a site (its own change there); one that carries the
// species' base blocks the site, and so do two dissenters. Below six the nine-in-ten rule needs all of them.
TEST(AncestrySites, ConsensusToleratesOneCongenerWithABaseOfItsOwn) {
    std::mt19937 rng(15);
    auto const own = RandomSequence(900, rng);
    auto make = [&](size_t n) {
        std::vector<std::string> seqs(n, own);
        for (auto& s : seqs) s[100] = Alt(own[100], 0);                        // all: a site
        for (size_t i = 1; i < n; i++) seqs[i][200] = Alt(own[200], 0);        // all but the first, which has the species'
        for (size_t i = 0; i < n; i++) seqs[i][300] = Alt(own[300], i == 0 ? 1 : 0);  // all but one, which has a third base
        for (size_t i = 2; i < n; i++) seqs[i][400] = Alt(own[400], i == 2 ? 1 : 0);  // one third base and two the species'
        for (size_t i = 2; i < n; i++) seqs[i][500] = Alt(own[500], i < 4 ? 1 : 0);   // two third bases, two the species'
        std::vector<an::Comparison> comparisons;
        for (size_t i = 0; i < n; i++) {
            auto comparison = an::CompareCopies(own, seqs[i]);
            EXPECT_TRUE(comparison.has_value());
            comparison->sites.congener = static_cast<uint32_t>(31 + i);
            comparisons.push_back(std::move(*comparison));
        }
        return an::Consensus(own.size(), comparisons);
    };
    auto const six = make(6);
    EXPECT_EQ(six.positions, (std::vector<uint16_t>{ 100, 300 }));
    EXPECT_EQ(six.bases[1], an::Code(Alt(own[300], 0)));  // the five's base, not the dissenter's
    auto const five = make(5);
    EXPECT_EQ(five.positions, (std::vector<uint16_t>{ 100 }));
    EXPECT_TRUE(an::Agree(5, 1, 6));
    EXPECT_FALSE(an::Agree(5, 0, 6));   // the sixth carries the species' base
    EXPECT_FALSE(an::Agree(4, 2, 6));
    EXPECT_FALSE(an::Agree(4, 1, 5));
    EXPECT_TRUE(an::Agree(9, 0, 10));   // nine in ten, as before
    EXPECT_TRUE(an::Agree(10, 1, 11));
}

// The indels alike: an indel site where nine in ten of the congeners compared at its flanks carry it (all of four), the
// nearest congener's with fewer than three compared.
TEST(AncestrySites, ConsensusTakesTheIndelsTheCongenersShare) {
    std::mt19937 rng(16);
    auto const own = RandomSequence(900, rng);
    std::vector<an::Comparison> comparisons;
    for (uint32_t taxid = 21; taxid <= 24; taxid++) {
        // Each differs at a base of its own, the fourth at two (the farthest: an indel does not count in the identity).
        auto other = Substituted(own, taxid == 24 ? std::vector<size_t>{ 250, 275 } : std::vector<size_t>{ 100 + 50 * (taxid - 21) }, rng);
        other.erase(400, 6);                                 // all four lack the species' bases 400-405
        if (taxid == 24) other.insert(600, "GATTACA");       // one of four has 7 bases the species lacks
        auto comparison = an::CompareCopies(own, other);
        ASSERT_TRUE(comparison.has_value());
        ASSERT_EQ(comparison->sites.indels.size(), taxid == 24 ? 2u : 1u);
        comparison->sites.congener = taxid;
        comparisons.push_back(std::move(*comparison));
    }
    auto const sites = an::Consensus(own.size(), comparisons);
    ASSERT_EQ(sites.indels.size(), 1u);
    EXPECT_EQ(sites.indels[0].length, 6);
    EXPECT_GE(sites.indels[0].position, 395u);
    EXPECT_LE(sites.indels[0].position, 405u);
    EXPECT_TRUE(sites.positions.empty());  // each base substitution is one congener's own
    // Two congeners: the nearest's indels, both of them for the one with the insertion if it is the nearest.
    auto const two = an::Consensus(own.size(), { comparisons[3], comparisons[0] });
    EXPECT_EQ(two.congener, 21u);  // one substitution, the highest identity
    EXPECT_EQ(two.indels.size(), 1u);
    auto const alone = an::Consensus(own.size(), { comparisons[3] });
    EXPECT_EQ(alone.indels.size(), 2u);
}

// A record counts an indel site it is aligned five bases beyond on both sides: the congeners' state when it has a gap
// of the site's kind and length within four bases of it, the species' when it has no gap near, neither otherwise.
TEST(AncestrySites, CountFindsTheIndelSites) {
    an::Sites sites;
    sites.indels = { { 40, 6 }, { 70, -5 } };  // the species' extra bases 40-45; 5 bases the congeners have before 70
    std::string const seq(100, 'A');
    auto c = an::Count(sites, "100M", 1, seq);
    EXPECT_EQ(c.indel_sites, 2u);
    EXPECT_EQ(c.indel_agree, 2u);
    EXPECT_EQ(c.indel_congener, 0u);
    EXPECT_EQ(c.sites, 0u);  // no base site
    c = an::Count(sites, "40M6D54M", 1, std::string(94, 'A'));  // deletes the species' extra bases
    EXPECT_EQ(c.indel_sites, 2u);
    EXPECT_EQ(c.indel_congener, 1u);
    EXPECT_EQ(c.indel_agree, 1u);
    c = an::Count(sites, "70M5I30M", 1, std::string(105, 'A'));  // inserts the congeners'
    EXPECT_EQ(c.indel_sites, 2u);
    EXPECT_EQ(c.indel_congener, 1u);
    c = an::Count(sites, "42M6D52M", 1, std::string(94, 'A'));  // the gap slid by two bases: still the congeners'
    EXPECT_EQ(c.indel_congener, 1u);
    c = an::Count(sites, "40M2D58M", 1, std::string(98, 'A'));  // another gap at the site: neither way
    EXPECT_EQ(c.indel_sites, 1u);
    EXPECT_EQ(c.indel_agree, 1u);
    c = an::Count(sites, "60M", 1, std::string(60, 'A'));  // aligned through the first site only
    EXPECT_EQ(c.indel_sites, 1u);
    c = an::Count(sites, "60M", 38, std::string(60, 'A'));  // starts too near the first site: the second only
    EXPECT_EQ(c.indel_sites, 1u);
    EXPECT_EQ(an::Count(sites, "100M", 1, "*").indel_sites, 2u);  // the indels need no sequence
    EXPECT_EQ(an::Count(an::Sites{}, "100M", 1, seq).indel_sites, 0u);
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

// Two congeners with the gene: fewer than three, so the sites are the nearest pairing one's differences.
TEST(AncestrySites, CacheTakesTheNearestCongenerWithThePairingCopy) {
    std::mt19937 rng(13);
    auto const own = RandomSequence(800, rng);
    FakeGenomes genomes;
    genomes.Set(1, 5, own);
    genomes.Set(3, 5, Substituted(own, { 50, 150, 250 }, rng));  // the second-nearest congener: the nearest with the gene
    genomes.Set(4, 5, Substituted(own, { 60, 160, 260, 360 }, rng));
    genomes.Set(1, 6, own);
    genomes.Set(2, 6, RandomSequence(800, rng));  // the nearest congener's copy of gene 6 does not pair
    genomes.Set(3, 6, Substituted(own, { 70 }, rng));
    sn::Table neighbours;
    neighbours.Set(1, { { 2, 0.02f }, { 3, 0.04f }, { 4, 0.06f } });
    an::Cache cache;
    auto const s5 = cache.Get(1, 5, genomes, neighbours);
    EXPECT_EQ(s5->congener, 3u);
    EXPECT_EQ(s5->congeners, 2u);
    EXPECT_EQ(s5->positions, (std::vector<uint16_t>{ 50, 150, 250 }));
    EXPECT_EQ(cache.Get(1, 5, genomes, neighbours).get(), s5.get());  // computed once
    auto const s6 = cache.Get(1, 6, genomes, neighbours);
    EXPECT_EQ(s6->congener, 3u);
    EXPECT_EQ(s6->congeners, 1u);
    EXPECT_EQ(s6->positions.size(), 1u);
    // A species without neighbours, or without the gene: empty; without the table: empty and nothing cached.
    EXPECT_TRUE(cache.Get(3, 5, genomes, neighbours)->Empty());
    EXPECT_TRUE(cache.Get(1, 7, genomes, neighbours)->Empty());
    EXPECT_EQ(cache.Size(), 4u);
    EXPECT_TRUE(cache.Get(1, 5, genomes, sn::Table{})->Empty());
    cache.Clear();
    EXPECT_EQ(cache.Size(), 0u);
}

// With three congeners or more the cache takes the consensus: the sites every one of them shares (nine in ten).
TEST(AncestrySites, CacheTakesTheConsensusOfTheCongeners) {
    std::mt19937 rng(15);
    auto const own = RandomSequence(800, rng);
    FakeGenomes genomes;
    genomes.Set(1, 5, own);
    std::string shared = own;
    shared[100] = Alt(own[100], 0);
    shared[500] = Alt(own[500], 1);
    genomes.Set(2, 5, Substituted(shared, { 200 }, rng));   // each differs at 100 and 500 alike, and at a site of its own
    genomes.Set(3, 5, Substituted(shared, { 300 }, rng));
    genomes.Set(4, 5, Substituted(shared, { 400 }, rng));
    genomes.Set(5, 5, RandomSequence(800, rng));            // does not pair: not a vote
    sn::Table neighbours;
    neighbours.Set(1, { { 2, 0.02f }, { 3, 0.03f }, { 4, 0.04f }, { 5, 0.05f } });
    an::Cache cache;
    auto const sites = cache.Get(1, 5, genomes, neighbours);
    EXPECT_EQ(sites->positions, (std::vector<uint16_t>{ 100, 500 }));
    EXPECT_EQ(sites->bases, (std::vector<uint8_t>{ an::Code(shared[100]), an::Code(shared[500]) }));
    EXPECT_EQ(sites->congeners, 3u);
    EXPECT_EQ(sites->congener, 2u);  // all three at the same identity: the first
}

// With congener_gaps.tsv the cache compares a copy with the gene's nearest congener by alignment too: also one that
// species_neighbours.tsv does not list (farther over all genes, or beyond its 16 within 0.15), and even without that
// table; without the gene's nearest (an older table, or none for the gene) the species' neighbours as before.
TEST(AncestrySites, CacheTakesTheGenesNearestCongenerFromTheGaps) {
    std::mt19937 rng(17);
    auto const own = RandomSequence(800, rng);
    FakeGenomes genomes;
    genomes.Set(1, 5, own);
    genomes.Set(2, 5, Substituted(own, { 40, 140, 240, 340, 440 }, rng));  // the species' nearest over all genes
    genomes.Set(7, 5, Substituted(own, { 100 }, rng));                     // gene 5's nearest, not a listed neighbour
    genomes.Set(1, 6, own);
    genomes.Set(2, 6, Substituted(own, { 30, 130 }, rng));
    sn::Table neighbours;
    neighbours.Set(1, { { 2, 0.02f } });
    auto const gaps = protal::congener_gaps::Table::FromRows({ { 1, 5, protal::congener_gaps::Gap{ 13, 63, 2, 7 } },
                                                               { 1, 6, protal::congener_gaps::Gap{ 25, 25, 1, 0 } } });
    an::Cache cache;
    auto const s5 = cache.Get(1, 5, genomes, neighbours, &gaps);
    EXPECT_EQ(s5->congener, 7u);  // the nearest of the two compared
    EXPECT_EQ(s5->congeners, 2u);
    EXPECT_EQ(s5->positions, std::vector<uint16_t>{ 100 });
    // Gene 6's entry names no congener: the species' neighbour.
    EXPECT_EQ(cache.Get(1, 6, genomes, neighbours, &gaps)->congener, 2u);
    // Without species_neighbours.tsv, the gaps' congener still serves.
    an::Cache alone;
    EXPECT_EQ(alone.Get(1, 5, genomes, sn::Table{}, &gaps)->congener, 7u);
    EXPECT_TRUE(alone.Get(1, 6, genomes, sn::Table{}, &gaps)->Empty());
    // Without the gaps: the species' nearest, as before.
    an::Cache before;
    EXPECT_EQ(before.Get(1, 5, genomes, neighbours)->congener, 2u);
    EXPECT_EQ(before.Get(1, 5, genomes, neighbours)->positions.size(), 5u);
}
