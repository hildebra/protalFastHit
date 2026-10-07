// Unit tests for what a taxon's call takes from the other taxa of its sample (SampleContext.h,
// MicrobialProfile::ApplySampleContext): sketched distances between references, the abundance-weighted assignment of
// ambiguous reads, the prior adjusted to a sample and the calls at a target share of false calls, the model header that
// holds them, and the relatives features and the singleton rule on a profile.
#include <gtest/gtest.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>
#include <unistd.h>
#include "Profiling/Profiler.h"
#include "Profiling/SampleContext.h"

using namespace protal;
namespace ctx = protal::profiler::context;

namespace {
    std::string RandomSequence(size_t length, std::mt19937& rng) {
        static constexpr char kBases[] = "ACGT";
        std::string seq(length, 'A');
        for (auto& c : seq) c = kBases[rng() % 4];
        return seq;
    }

    // seq with a share `rate` of its bases changed to another base.
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

    SamEntry MakeSam(uint32_t taxid, uint32_t geneid, std::string seq, POS_t pos, std::string alternatives = "*") {
        SamEntry sam;
        sam.m_qname = "read";
        sam.m_rname = std::to_string(taxid) + "_" + std::to_string(geneid);
        sam.m_pos = pos;
        sam.m_mapq = 60;
        sam.m_cigar = std::to_string(seq.size()) + "M";
        sam.m_qual = std::string(seq.size(), 'I');
        sam.m_seq = std::move(seq);
        sam.m_alternatives = std::move(alternatives);
        return sam;
    }

    // A reference of taxa with genes 1..n each: reference.fna and reference.map in a folder of their own.
    struct Reference {
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;
        std::map<uint32_t, std::vector<std::string>> genes;  // taxid -> gene i at [i - 1]

        explicit Reference(std::map<uint32_t, std::vector<std::string>> taxa) : genes(std::move(taxa)) {
            dir = std::filesystem::temp_directory_path() / ("protal_sample_context_" + std::to_string(::getpid()) + "_" +
                                                            std::to_string(reinterpret_cast<uintptr_t>(this)));
            std::filesystem::create_directories(dir);
            std::ofstream fna(dir / "reference.fna");
            std::ofstream map(dir / "reference.map");
            size_t offset = 0;
            for (auto const& [taxid, seqs] : genes) {
                for (size_t i = 0; i < seqs.size(); i++) {
                    std::string const header = ">" + std::to_string(taxid) + "_" + std::to_string(i + 1) + "\n";
                    fna << header << seqs[i] << '\n';
                    map << taxid << '\t' << i + 1 << '\t' << offset + header.size() << '\t'
                        << offset + header.size() + seqs[i].size() << '\n';
                    offset += header.size() + seqs[i].size() + 1;
                }
            }
            fna.close();
            map.close();
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
            loader->LoadAllGenomes();
        }
        ~Reference() { std::filesystem::remove_all(dir); }
    };
}

namespace {
    // GeneSketch as first written: every distinct k-mer collected and sorted, hashed, the smallest kept.
    std::vector<uint32_t> ReferenceGeneSketch(std::string_view seq) {
        auto const kmers = gene_conservation::Kmers(seq);
        std::vector<uint32_t> hashes;
        hashes.reserve(kmers.size());
        for (uint32_t const kmer : kmers) hashes.push_back(ctx::MixKmer(kmer));
        size_t const keep = std::min(ctx::kSketchSize, hashes.size());
        std::partial_sort(hashes.begin(), hashes.begin() + static_cast<std::ptrdiff_t>(keep), hashes.end());
        hashes.resize(keep);
        return hashes;
    }
}

// Random genes, genes with Ns, short ones, repeats of a few k-mers (fewer distinct than the sketch holds, which sorts
// everything), half-repeated genes and genes of one base: the same hashes as the reference, in the same order.
TEST(SampleContext, SketchesEqualTheReferencesToTheLastHash) {
    std::mt19937 rng(7);
    std::vector<std::string> genes;
    for (int i = 0; i < 30; i++) genes.push_back(RandomSequence(200 + rng() % 3000, rng));
    for (int i = 0; i < 10; i++) {
        auto gene = RandomSequence(1000, rng);
        for (int n = 0; n < 20; n++) gene[rng() % gene.size()] = 'N';
        genes.push_back(gene);
    }
    for (size_t length : { 0, 5, 11, 12, 13, 40, 75, 76, 100 }) genes.push_back(RandomSequence(length, rng));
    std::string const unit = RandomSequence(30, rng);
    std::string repeated;
    for (int i = 0; i < 40; i++) repeated += unit;  // 30 distinct k-mers, 1,200 in all
    genes.push_back(repeated);
    genes.push_back(repeated.substr(0, 600) + RandomSequence(600, rng));
    genes.push_back(std::string(500, 'A'));
    genes.push_back(std::string(500, 'A') + RandomSequence(300, rng));
    genes.push_back("acgtacgtnnACGT" + RandomSequence(200, rng));
    for (auto const& gene : genes) {
        auto const expected = ReferenceGeneSketch(gene);
        auto const got = ctx::GeneSketch(gene);
        EXPECT_EQ(got, expected) << "gene of " << gene.size() << " bases";
    }
    EXPECT_LT(ctx::GeneSketch(repeated).size(), ctx::kSketchSize);
}

TEST(SampleContext, SketchesOfIdenticalGenesAreAtDistanceZeroAndOfUnrelatedOnesFar) {
    std::mt19937 rng(1);
    auto const gene = RandomSequence(1000, rng);
    auto const sketch = ctx::GeneSketch(gene);
    ASSERT_EQ(sketch.size(), ctx::kSketchSize);
    EXPECT_TRUE(std::is_sorted(sketch.begin(), sketch.end()));
    EXPECT_EQ(ctx::SketchDistance(sketch, sketch), 0.0);
    EXPECT_EQ(ctx::SketchDistance(sketch, ctx::GeneSketch(RandomSequence(1000, rng))), 1.0);
    EXPECT_EQ(ctx::SketchDistance(sketch, {}), 1.0);
    // 5% of the bases changed: about Mash's distance on all k-mers.
    auto const mutated = Mutated(gene, 0.05, rng);
    double const exact = gene_conservation::MashDistance(gene_conservation::Kmers(gene), gene_conservation::Kmers(mutated));
    double const sketched = ctx::SketchDistance(sketch, ctx::GeneSketch(mutated));
    EXPECT_NEAR(sketched, exact, 0.025);
    EXPECT_GT(sketched, 0.02);
}

TEST(SampleContext, TwoReferencesNeedTenSharedGenesForADistance) {
    std::mt19937 rng(2);
    ctx::TaxonSketch a, b;
    for (uint32_t g = 1; g <= 9; g++) {
        auto const gene = RandomSequence(300, rng);
        a.emplace_back(g, ctx::GeneSketch(gene));
        b.emplace_back(g, ctx::GeneSketch(gene));
    }
    auto nine = ctx::SketchedTaxonDistances(a, b);
    EXPECT_EQ(nine.distance, ctx::kFarDistance);
    EXPECT_EQ(nine.genes.size(), 9u);
    auto const gene = RandomSequence(300, rng);
    a.emplace_back(10, ctx::GeneSketch(gene));
    b.emplace_back(10, ctx::GeneSketch(gene));
    b.emplace_back(11, ctx::GeneSketch(RandomSequence(300, rng)));  // only in b: not compared
    auto ten = ctx::SketchedTaxonDistances(a, b);
    EXPECT_EQ(ten.distance, 0.0);
    EXPECT_EQ(ten.genes.size(), 10u);
}

TEST(SampleContext, SpillRateFallsTenfoldPerDecade) {
    EXPECT_DOUBLE_EQ(ctx::SpillRate(0), ctx::kSpillAtZero);
    EXPECT_NEAR(ctx::SpillRate(ctx::kSpillDecade), ctx::kSpillAtZero / 10, 1e-15);
    EXPECT_NEAR(ctx::SpillRate(2 * ctx::kSpillDecade), ctx::kSpillAtZero / 100, 1e-15);
}

TEST(SampleContext, TheSamplesComplexityIsItsTaxaTheirLowIdentityShareAndTheirMedianIdentity) {
    // Taxa (fragments, identity, low-identity share) in taxid order; one without fragments counts for nothing.
    std::vector<ctx::TaxonSummary> taxa = { { 30, 0.99, 0.0 }, { 0, 0.5, 1.0 }, { 10, 0.95, 0.5 }, { 60, 0.97, 0.1 },
                                            { 2, 0.80, 1.0 } };
    auto const c = ctx::SampleComplexityOf(taxa);
    EXPECT_NEAR(c.log_taxa, std::log10(4.0), 1e-12);
    EXPECT_NEAR(c.low_identity, (10 * 0.5 + 60 * 0.1 + 2 * 1.0) / 102, 1e-12);
    // The median of the taxa with at least 10 fragments (100 of them): 10 at 0.95, then 60 at 0.97 reach half.
    EXPECT_EQ(c.identity, 0.97);
    // With half reached exactly at a taxon, that taxon's identity (numpy's searchsorted of the cumulative sum).
    EXPECT_EQ(ctx::SampleComplexityOf({ { 20, 0.9, 0 }, { 20, 0.99, 0 } }).identity, 0.9);
    // The taxa are sorted by identity whatever order they come in.
    EXPECT_EQ(ctx::SampleComplexityOf({ { 40, 0.99, 0 }, { 15, 0.90, 0 }, { 30, 0.95, 0 } }).identity, 0.95);
    // No taxon of 10 fragments: the median over every taxon with fragments, here the one holding 3 of the 5.
    auto const thin = ctx::SampleComplexityOf({ { 1, 0.98, 0 }, { 3, 0.93, 0.25 }, { 1, 1.0, 0 } });
    EXPECT_EQ(thin.identity, 0.93);
    EXPECT_NEAR(thin.low_identity, 0.75 / 5, 1e-12);
    // An empty sample: zeros.
    auto const empty = ctx::SampleComplexityOf({ { 0, 0.9, 0.5 } });
    EXPECT_EQ(empty.log_taxa, 0.0);
    EXPECT_EQ(empty.low_identity, 0.0);
    EXPECT_EQ(empty.identity, 0.0);
}

TEST(SampleContext, TheEditRatioIsTheTaxonsOwnDivergence) {
    EXPECT_EQ(ctx::EditRatio({ 10, 10, 0, 0 }), ctx::kEditRatio);  // no aligned bases
    EXPECT_NEAR(ctx::EditRatio({ 10, 10, 2, 100 }), 0.02 / 0.98, 1e-12);
    EXPECT_NEAR(ctx::EditRatio({ 10, 10, 0, 100 }), ctx::kMinDivergence / (1 - ctx::kMinDivergence), 1e-12);
    EXPECT_NEAR(ctx::EditRatio({ 10, 10, 90, 100 }), ctx::kMaxDivergence / (1 - ctx::kMaxDivergence), 1e-12);
}

// A taxon whose reads all fit a congener a hundred times as abundant one edit worse holds the congener's reads; one
// whose reads fit it five edits worse, or have no alternative, keeps them.
TEST(SampleContext, AnAbundantCongenerTakesTheReadsItFitsNearlyAsWell) {
    std::map<uint32_t, ctx::RecordCounts> counts = {
        { 1, { 1000, 1000, 5000, 150000 } },  // 3.3% divergence: edit ratio 0.034
        { 2, { 10, 10, 50, 1500 } },
        { 3, { 10, 10, 50, 1500 } },
        { 4, { 10, 6, 50, 1500 } } };
    ctx::AmbiguityClasses classes;
    classes[{ 2, 1, 1, 1 }] = 10;  // every record of 2 fits 1 one edit worse
    classes[{ 4, 1, 1, 5 }] = 6;   // the kept records of 4 fit 1 five edits worse
    auto const shares = ctx::AbundanceWeightedShares(classes, counts);
    EXPECT_LT(shares.at(2).all, 0.35);
    EXPECT_EQ(shares.at(2).all, shares.at(2).kept);
    EXPECT_EQ(shares.at(3).all, 1.0);
    EXPECT_GT(shares.at(4).kept, 0.99);
    EXPECT_EQ(shares.at(1).all, 1.0);  // its reads have no alternatives
    // The same counts give the same shares, whatever order the classes were noted in.
    auto const again = ctx::AbundanceWeightedShares(classes, counts);
    EXPECT_EQ(again.at(2).all, shares.at(2).all);

    // A tie (0 edits more) goes almost all to the abundant taxon.
    ctx::AmbiguityClasses ties;
    ties[{ 2, 1, 1, 0 }] = 10;
    EXPECT_LT(ctx::AbundanceWeightedShares(ties, counts).at(2).all, 0.05);
    // With few reads of its own the congener takes less: two taxa of 10 records share their ties evenly.
    std::map<uint32_t, ctx::RecordCounts> even = { { 1, { 10, 10, 50, 1500 } }, { 2, { 10, 10, 50, 1500 } } };
    ctx::AmbiguityClasses both;
    both[{ 1, 1, 2, 0 }] = 10;
    both[{ 2, 1, 1, 0 }] = 10;
    auto const split = ctx::AbundanceWeightedShares(both, even);
    EXPECT_NEAR(split.at(1).all, 0.5, 1e-9);
    EXPECT_NEAR(split.at(2).all, 0.5, 1e-9);
}

namespace {
    // AbundanceWeightedShares as first written (be35d15): hash maps and a pow per visit. The version in
    // SampleContext.h runs on indices and must give the same shares to the last bit.
    std::unordered_map<uint32_t, ctx::OwnShares> ReferenceShares(ctx::AmbiguityClasses const& classes,
                                                                 std::map<uint32_t, ctx::RecordCounts> const& counts,
                                                                 double tolerance = 1e-10) {
        std::unordered_map<uint32_t, double> plain, weight, ratio;
        for (auto const& [taxid, c] : counts) {
            plain[taxid] = static_cast<double>(c.records);
            weight[taxid] = static_cast<double>(c.records);
            ratio[taxid] = ctx::EditRatio(c);
        }
        for (auto const& [key, n] : classes) plain[key[0]] -= static_cast<double>(n);
        std::vector<double> shares;
        auto posterior = [&](ctx::AmbiguityKey const& key) {
            shares.assign(1 + (key.size() - 2) / 2, 0.0);
            auto const own = weight.find(key[0]);
            shares[0] = own == weight.end() ? 0 : own->second;
            double total = shares[0];
            for (size_t i = 2, a = 1; i + 1 < key.size(); i += 2, a++) {
                auto const found = weight.find(key[i]);
                shares[a] = found == weight.end() ? 0 : found->second * std::pow(ratio.at(key[i]), static_cast<double>(key[i + 1]));
                total += shares[a];
            }
            if (total <= 0) {
                shares.assign(shares.size(), 0.0);
                shares[0] = 1;
                return;
            }
            for (auto& s : shares) s /= total;
        };
        for (size_t iteration = 0; iteration < ctx::kEmIterations && !classes.empty(); iteration++) {
            std::unordered_map<uint32_t, double> next = plain;
            for (auto const& [key, n] : classes) {
                posterior(key);
                next[key[0]] += static_cast<double>(n) * shares[0];
                for (size_t i = 2, a = 1; i + 1 < key.size(); i += 2, a++) {
                    if (shares[a] > 0) next[key[i]] += static_cast<double>(n) * shares[a];
                }
            }
            double change = 0;
            for (auto const& [taxid, w] : next) {
                auto const before = weight.find(taxid);
                double const old = before == weight.end() ? 0 : before->second;
                change = std::max(change, std::abs(w - old) / std::max(1.0, old));
            }
            weight = std::move(next);
            if (change < tolerance) break;
        }
        std::unordered_map<uint32_t, double> own_all, own_kept;
        for (auto const& [key, n] : classes) {
            posterior(key);
            own_all[key[0]] += static_cast<double>(n) * shares[0];
            if (key[1]) own_kept[key[0]] += static_cast<double>(n) * shares[0];
        }
        std::unordered_map<uint32_t, uint64_t> ambiguous_kept;
        for (auto const& [key, n] : classes) {
            if (key[1]) ambiguous_kept[key[0]] += n;
        }
        std::unordered_map<uint32_t, ctx::OwnShares> result;
        for (auto const& [taxid, c] : counts) {
            ctx::OwnShares s;
            if (c.records > 0) s.all = std::clamp((plain[taxid] + own_all[taxid]) / static_cast<double>(c.records), 0.0, 1.0);
            if (c.kept > 0) {
                double const plain_kept = static_cast<double>(c.kept) - static_cast<double>(ambiguous_kept[taxid]);
                s.kept = std::clamp((plain_kept + own_kept[taxid]) / static_cast<double>(c.kept), 0.0, 1.0);
            }
            result[taxid] = s;
        }
        return result;
    }
}

// Random samples: taxa of very different depths and divergences, classes with 0-4 alternatives at 0-5 edits, ties,
// alternatives the counts lack (a ZA taxon without a best record) and an owner the counts lack; the shares of every
// taxon equal the reference's exactly (EXPECT_EQ on the doubles), on every sample. Both run every sweep (tolerance
// 0), as the two stop by different rules (StopsWhenTheSharesAreStable).
TEST(SampleContext, SharesOnIndicesEqualTheReferenceToTheLastBit) {
    std::mt19937 rng(20261003);
    size_t compared = 0;
    for (int sample = 0; sample < 40; sample++) {
        size_t const taxa = 2 + rng() % 60;
        std::map<uint32_t, ctx::RecordCounts> counts;
        for (size_t t = 0; t < taxa; t++) {
            uint32_t const taxid = 1 + static_cast<uint32_t>(rng() % 200);
            uint64_t const records = 1 + (rng() % 7 == 0 ? rng() % 20000 : rng() % 50);
            uint64_t const kept = rng() % (records + 1);
            uint64_t const aligned = records * 150;
            uint64_t const differences = static_cast<uint64_t>(static_cast<double>(aligned) * (rng() % 1000) / 1000.0 * 0.1);
            counts[taxid] = { records, kept, differences, rng() % 5 == 0 ? 0 : aligned };
        }
        std::vector<uint32_t> ids;
        for (auto const& [taxid, _] : counts) ids.push_back(taxid);
        ctx::AmbiguityClasses classes;
        size_t const n_classes = rng() % 300;
        for (size_t k = 0; k < n_classes; k++) {
            ctx::AmbiguityKey key;
            key.push_back(rng() % 50 == 0 ? 999 : ids[rng() % ids.size()]);  // now and then an owner without records
            key.push_back(rng() % 2);
            size_t const alternatives = rng() % 5;
            for (size_t a = 0; a < alternatives; a++) {
                key.push_back(rng() % 20 == 0 ? 998 : ids[rng() % ids.size()]);  // or an alternative without records
                key.push_back(rng() % 6);
            }
            classes[key] += 1 + rng() % 30;
        }
        auto const expected = ReferenceShares(classes, counts, 0);
        auto const got = ctx::AbundanceWeightedShares(classes, counts, 0);
        ASSERT_EQ(got.size(), expected.size());
        for (auto const& [taxid, s] : expected) {
            ASSERT_TRUE(got.contains(taxid));
            EXPECT_EQ(got.at(taxid).all, s.all) << "sample " << sample << " taxon " << taxid;
            EXPECT_EQ(got.at(taxid).kept, s.kept) << "sample " << sample << " taxon " << taxid;
            compared++;
        }
    }
    EXPECT_GT(compared, 500u);
}

// The sweeps stop once no taxon's own share moved by kEmTolerance in a sweep: a sample that converges stops early
// with the shares of all kEmIterations sweeps to that tolerance; a taxon that loses its last reads to a taxon which
// explains them exactly as well (an abundant congener whose weight times its edit factor equals the taxon's reads)
// loses them ever more slowly and runs to the cap, where both stop alike.
TEST(SampleContext, StopsWhenTheSharesAreStable) {
    std::map<uint32_t, ctx::RecordCounts> counts = {
        { 1, { 1000, 1000, 5000, 150000 } }, { 2, { 50, 50, 250, 7500 } }, { 3, { 30, 30, 150, 4500 } }, { 4, { 20, 20, 100, 3000 } } };
    ctx::AmbiguityClasses classes;
    classes[{ 2, 1, 1, 1 }] = 40;
    classes[{ 2, 1, 3, 0 }] = 10;
    classes[{ 3, 1, 2, 1, 1, 2 }] = 20;
    classes[{ 4, 1, 1, 2, 3, 0 }] = 15;
    classes[{ 1, 1, 2, 3 }] = 100;
    auto const stopped = ctx::AbundanceWeightedShares(classes, counts);
    auto const full = ctx::AbundanceWeightedShares(classes, counts, 0);
    auto const few = ctx::AbundanceWeightedShares(classes, counts, 0, 2);
    for (auto const& [taxid, s] : full) {
        EXPECT_NEAR(stopped.at(taxid).all, s.all, 1e-5) << "taxon " << taxid;
        EXPECT_NEAR(stopped.at(taxid).kept, s.kept, 1e-5) << "taxon " << taxid;
    }
    // Two sweeps are not enough on this sample: the rule did not stop that early.
    bool moved = false;
    for (auto const& [taxid, s] : full) moved |= std::abs(few.at(taxid).all - s.all) > 1e-5;
    EXPECT_TRUE(moved);

    // A tie: taxon 6's 10 reads fit taxon 5 two edits worse, and taxon 5's 4,000 records at a factor of 0.05^2 are
    // 10 effective reads, as many as taxon 6 has. Taxon 6's weight w goes to 10 w / (w + 10): harmonically
    // (10 / (k + 1) after k sweeps) while it is near 10, then at a rate within 0.2% of standing still (taxon 5
    // gained up to 5 reads); after kEmIterations sweeps its share is still about a hundredth and moving by more than
    // the tolerance: the rule runs to the cap, where the shares equal the full run's exactly.
    std::map<uint32_t, ctx::RecordCounts> tie = { { 5, { 4000, 4000, 0, 0 } }, { 6, { 10, 10, 0, 0 } } };  // ratio kEditRatio
    ctx::AmbiguityClasses drain;
    drain[{ 6, 1, 5, 2 }] = 10;
    auto const capped = ctx::AbundanceWeightedShares(drain, tie);
    auto const capped_full = ctx::AbundanceWeightedShares(drain, tie, 0);
    EXPECT_EQ(capped.at(6).all, capped_full.at(6).all);  // both ran every sweep
    EXPECT_GT(capped.at(6).all, 1e-3);
    EXPECT_LT(capped.at(6).all, 0.05);
    auto const capped_half = ctx::AbundanceWeightedShares(drain, tie, 0, ctx::kEmIterations / 2);
    EXPECT_GT(std::abs(capped_half.at(6).all - capped.at(6).all), 1e-3);  // still moving at half the cap
}

TEST(SampleContext, TheSamplesPriorFollowsItsCandidates) {
    // Probabilities whose mean is the prior are kept.
    std::vector<double> const q = { 0.9, 0.1, 0.5, 0.5 };
    double rate = 0;
    auto const same = ctx::SampleAdjusted(q, 0.5, &rate);
    EXPECT_NEAR(rate, 0.5, 1e-9);
    for (size_t i = 0; i < q.size(); i++) EXPECT_NEAR(same[i], q[i], 1e-9);
    // Many unlikely candidates lower the sample's rate and every probability with it.
    std::vector<double> deep = { 0.95, 0.9 };
    deep.resize(200, 0.05);
    auto const adjusted = ctx::SampleAdjusted(deep, 0.3, &rate);
    EXPECT_LT(rate, 0.1);
    EXPECT_LT(adjusted[0], 0.95);
    EXPECT_LT(adjusted[5], 0.05);
}

TEST(SampleContext, CallsKeepTheExpectedShareOfFalseOnesAtTheTarget) {
    ctx::CalibrationCurve const identity = { { 0, 0 }, { 1, 1 } };
    EXPECT_EQ(ctx::Calibrated(identity, 0.3), 0.3);
    EXPECT_EQ(ctx::Calibrated({ { 0.2, 0.1 }, { 0.8, 0.7 } }, 0.1), 0.1);  // the end points' beyond them
    EXPECT_EQ(ctx::Calibrated({ { 0.2, 0.1 }, { 0.8, 0.7 } }, 0.9), 0.7);
    EXPECT_NEAR(ctx::Calibrated({ { 0.2, 0.1 }, { 0.8, 0.7 } }, 0.5), 0.4, 1e-12);

    std::vector<double> const scores = { 0.99, 0.97, 0.9, 0.6, 0.3, 0.1, 0.05 };
    size_t last = 0;
    for (double fdr : { 0.001, 0.01, 0.05, 0.1, 0.2, 0.4 }) {
        auto const calls = ctx::FalseCallKnob(scores, identity, 0.5, fdr);
        EXPECT_GE(calls.called, last) << fdr;
        last = calls.called;
        EXPECT_LE(calls.expected_false, fdr * static_cast<double>(calls.called) + 1e-12);
        if (calls.called > 0) EXPECT_EQ(calls.knob, scores[calls.called - 1]);
    }
    EXPECT_GE(last, 4u);
    auto const none = ctx::FalseCallKnob(scores, identity, 0.5, 1e-9);
    EXPECT_EQ(none.called, 0u);
    EXPECT_GT(none.knob, 1.0);
    auto const empty = ctx::FalseCallKnob({}, identity, 0.5, 0.1);
    EXPECT_EQ(empty.called, 0u);
    // Tied scores are called together (every taxon at the knob's score).
    auto const tied = ctx::FalseCallKnob({ 0.95, 0.95, 0.1 }, identity, 0.5, 0.1);
    EXPECT_EQ(tied.knob, 0.95);
}

TEST(SampleContext, AModelsCalibratedCallsAreReadFromItsHeader) {
    auto const header = [](std::string extensions) {
        return "<PMML><Header description=\"x\">" + extensions + "<Application name=\"a\"/></Header><DataDictionary/></PMML>";
    };
    profiler::FalseCallModel model;
    EXPECT_EQ(profiler::ParseFalseCalls(header(""), model), "");
    EXPECT_TRUE(model.curve.empty());
    std::string const ok = "<Extension name=\"protal_calibration\" value=\"0:0.001,0.5:0.2,1:0.99\"/>"
                           "<Extension name=\"protal_prior\" value=\"0.21\"/><Extension name=\"protal_fdr\" value=\"0.05\"/>";
    ASSERT_EQ(profiler::ParseFalseCalls(header(ok), model), "");
    ASSERT_EQ(model.curve.size(), 3u);
    EXPECT_EQ(model.curve[1], (std::pair<double, double>{ 0.5, 0.2 }));
    EXPECT_EQ(model.prior, 0.21);
    EXPECT_EQ(model.fdr, 0.05);
    // Not all three, a decreasing probability, a score out of order or range, a prior or target outside (0, 1).
    for (std::string bad : { "<Extension name=\"protal_calibration\" value=\"0:0.1,1:0.9\"/>",
                             "<Extension name=\"protal_calibration\" value=\"0:0.5,1:0.4\"/><Extension name=\"protal_prior\" value=\"0.2\"/><Extension name=\"protal_fdr\" value=\"0.1\"/>",
                             "<Extension name=\"protal_calibration\" value=\"0.5:0.1,0.4:0.2\"/><Extension name=\"protal_prior\" value=\"0.2\"/><Extension name=\"protal_fdr\" value=\"0.1\"/>",
                             "<Extension name=\"protal_calibration\" value=\"0:0.1,1.5:0.2\"/><Extension name=\"protal_prior\" value=\"0.2\"/><Extension name=\"protal_fdr\" value=\"0.1\"/>",
                             "<Extension name=\"protal_calibration\" value=\"0:0.1,1:0.2\"/><Extension name=\"protal_prior\" value=\"1\"/><Extension name=\"protal_fdr\" value=\"0.1\"/>",
                             "<Extension name=\"protal_calibration\" value=\"0:0.1,1:0.2\"/><Extension name=\"protal_prior\" value=\"0.2\"/><Extension name=\"protal_fdr\" value=\"0\"/>",
                             "<Extension name=\"protal_calibration\" value=\"\"/><Extension name=\"protal_prior\" value=\"0.2\"/><Extension name=\"protal_fdr\" value=\"0.1\"/>" }) {
        EXPECT_NE(profiler::ParseFalseCalls(header(bad), model), "") << bad;
        EXPECT_TRUE(model.curve.empty()) << bad;
    }
}

// Taxa 1 and 2 are congeners (genus 10), their references 1% apart on genes 1-6 and 6% on genes 7-12; taxon 3 is of
// another genus of the family. Taxon 1 has 200 fragments.
TEST(SampleContext, AProfilesTaxaAreJudgedAgainstTheirAbundantRelatives) {
    std::mt19937 rng(5);
    std::vector<std::string> one, two, three;
    for (int g = 0; g < 12; g++) {
        one.push_back(RandomSequence(300, rng));
        two.push_back(Mutated(one.back(), g < 6 ? 0.01 : 0.06, rng));
        three.push_back(RandomSequence(300, rng));
    }
    Reference ref({ { 1, one }, { 2, two }, { 3, three } });
    auto genera = std::make_shared<std::vector<uint32_t>>(30, 0);
    auto families = std::make_shared<std::vector<uint32_t>>(30, 0);
    (*genera)[1] = (*genera)[2] = 10;
    (*genera)[3] = 11;
    (*families)[1] = (*families)[2] = (*families)[3] = 20;

    auto fill = [&](profiler::MicrobialProfile& profile, std::vector<uint32_t> const& two_genes, std::string const& za = "1:1") {
        profile.SetGenera(genera);
        profile.SetFamilies(families);
        int read = 0;
        for (int f = 0; f < 200; f++) {
            uint32_t const gene = static_cast<uint32_t>(f % 12) + 1;
            auto const sam = MakeSam(1, gene, one[gene - 1].substr(10, 100), 11);
            auto noted = sam;  // the abundant taxon's reads differ from its reference at 3% of their bases
            noted.m_cigar = "47M3X50M";
            profile.NoteRecord(1, gene, noted);
            ASSERT_TRUE(profile.AddSam(1, static_cast<int>(gene), sam, 1.0, true, read++));
        }
        for (auto gene : two_genes) {
            auto const sam = MakeSam(2, gene, two[gene - 1].substr(10, 100), 11, za);
            profile.NoteRecord(2, gene, sam);
            ASSERT_TRUE(profile.AddSam(2, static_cast<int>(gene), sam, 1.0, true, read++));
        }
        auto const sam = MakeSam(3, 1, three[0].substr(10, 100), 11);
        profile.NoteRecord(3, 1, sam);
        ASSERT_TRUE(profile.AddSam(3, 1, sam, 1.0, true, read++));
        profile.ApplyRecordEvidence();
    };

    profiler::MicrobialProfile profile(*ref.loader);
    profile.SetSingletonCongener(100);  // the rule is off by default
    fill(profile, { 1 });
    auto const& s1 = profile.GetTaxa().at(1).GetSampleEvidence();
    auto const& s2 = profile.GetTaxa().at(2).GetSampleEvidence();
    auto const& s3 = profile.GetTaxa().at(3).GetSampleEvidence();
    // Taxon 2: one fragment beside a congener of 200, which its read fits one edit worse.
    EXPECT_TRUE(profile.GetTaxa().at(2).Vetoed());
    EXPECT_EQ(s2.genus_top, 200u);
    EXPECT_NEAR(s2.genus_skew, std::log10(2.0 / 201), 1e-12);
    EXPECT_NEAR(s2.relative_skew, std::log10(2.0 / 201), 1e-12);
    EXPECT_NEAR(s2.genus_share, 1.0 / 201, 1e-12);
    EXPECT_GT(s2.relative_distance, 0.005);
    EXPECT_LT(s2.relative_distance, 0.08);
    EXPECT_NEAR(s2.relative_spill, std::log10(1.5 / (200 * ctx::SpillRate(s2.relative_distance) + 0.5)), 1e-12);
    EXPECT_LT(s2.em_own_share, 0.5);
    // Taxon 1 has no congener with more fragments; taxon 3 none at all, but one of 200 in another genus of its family.
    EXPECT_FALSE(profile.GetTaxa().at(1).Vetoed());
    EXPECT_NEAR(s1.genus_share, 200.0 / 201, 1e-12);
    EXPECT_EQ(s1.relative_distance, 1.0);
    EXPECT_EQ(s1.relative_close_share, 1.0);
    EXPECT_EQ(s1.em_own_share, 1.0);
    EXPECT_FALSE(profile.GetTaxa().at(3).Vetoed());
    EXPECT_EQ(s3.genus_top, 0u);
    EXPECT_NEAR(s3.genus_skew, std::log10(2.0), 1e-12);
    EXPECT_NEAR(s3.family_skew, std::log10(2.0 / 201), 1e-12);
    EXPECT_EQ(s3.relative_distance, 1.0);
    // The features carry them.
    std::map<std::string, double> features;
    for (auto const& [name, value] : profiler::TaxonFeatures(profile.GetTaxa().at(2))) features[name] = value;
    EXPECT_EQ(features.at("genus_top_fragments"), 200.0);
    EXPECT_EQ(features.at("genus_skew"), s2.genus_skew);
    EXPECT_EQ(features.at("relative_close_share"), s2.relative_close_share);

    // Without the singleton rule nothing is vetoed.
    profiler::MicrobialProfile no_rule(*ref.loader);
    no_rule.SetSingletonCongener(0);
    fill(no_rule, { 1 });
    EXPECT_FALSE(no_rule.GetTaxa().at(2).Vetoed());
    // Nor a single read that fits only its own reference (no alternative, identity 1): a minor congener's own read.
    profiler::MicrobialProfile own(*ref.loader);
    fill(own, { 1 }, "*");
    EXPECT_EQ(own.GetTaxa().at(2).GetSampleEvidence().em_own_share, 1.0);
    EXPECT_FALSE(own.GetTaxa().at(2).Vetoed());

    // Reads on the genes where the two references are most alike: about twice their share of the length; on the others none.
    profiler::MicrobialProfile close(*ref.loader);
    fill(close, { 1, 2, 3, 4 });
    EXPECT_FALSE(close.GetTaxa().at(2).Vetoed());
    EXPECT_NEAR(close.GetTaxa().at(2).GetSampleEvidence().relative_close_share, 2.0, 1e-9);
    profiler::MicrobialProfile far(*ref.loader);
    fill(far, { 8, 9, 10, 11 });
    EXPECT_NEAR(far.GetTaxa().at(2).GetSampleEvidence().relative_close_share, 0.0, 1e-9);
}
