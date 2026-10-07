// Unit tests for the per-copy tables and the candidates of 2026-10-07 (docs/claude/2026-10-07-congener-gaps): the gaps
// to the congeners' copies (CongenerGaps.h: the alignment distance, the scan of a gene's copies within their genera,
// congener_gaps.tsv), the foreign rates of a tiled scan (ForeignRatesTable.h), the untried candidates (ZC) and the
// adaptive candidates of the alignment handler, and the profiler's gaps, foreign and untried features.
#include <gtest/gtest.h>
#include <cmath>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "Core/AlignmentStrategy.h"
#include "IO/AlignmentOutputHandler.h"
#include "Profiling/Profiler.h"
#include "SequenceUtils/CongenerGaps.h"
#include "SequenceUtils/ForeignRatesTable.h"
#include "TestReference.h"

using namespace protal;
namespace cg = protal::congener_gaps;
namespace fr = protal::foreign_rates;

namespace {
    // seq with `n` substitutions at evenly spaced positions from `first` (every `step`-th base).
    std::string Substituted(std::string seq, size_t n, size_t first = 5, size_t step = 9) {
        for (size_t k = 0; k < n; k++) {
            size_t const i = first + k * step;
            seq[i] = seq[i] == 'A' ? 'C' : 'A';
        }
        return seq;
    }

    std::map<std::string, double> Features(profiler::Taxon const& taxon) {
        std::map<std::string, double> features;
        for (auto const& [name, value] : profiler::TaxonFeatures(taxon)) features[name] = value;
        return features;
    }

    // A record of the gene's own bases from 0-based `start`, MAPQ 60, with a CIGAR (default: all matches) and tags.
    std::string Record(test::LoadedReference const& ref, std::string const& qname, uint32_t taxid, uint32_t gene, int start,
                       int length, std::string const& cigar = "", std::string const& tags = "ZA:Z:*") {
        std::string const seq = ref.Gene(taxid, gene).substr(static_cast<size_t>(start), static_cast<size_t>(length));
        return qname + "\t0\t" + std::to_string(taxid) + "_" + std::to_string(gene) + '\t' + std::to_string(start + 1) + "\t60\t" +
               (cigar.empty() ? std::to_string(length) + "=" : cigar) + "\t*\t0\t0\t" + seq + '\t' + std::string(seq.size(), 'I') +
               '\t' + tags + '\n';
    }

    // Profiles a SAM (its text) as a run does, taxa 1 and 2 in genus 10, taxon 3 in genus 11.
    profiler::MicrobialProfile Profile(test::LoadedReference const& ref, std::string const& sam, size_t threads = 1) {
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
}

TEST(CongenerGaps, TheAlignedDistanceIsTheShareOfDifferingColumns) {
    std::mt19937 rng(3);
    WFA2Wrapper2 aligner(4, 6, 2, 0);
    std::string cigar;
    std::string const a = test::RandomSequence(900, rng);
    EXPECT_EQ(cg::AlignedDistance(a, a, aligner, cigar), 0.0);
    // 45 substitutions of 900 bases: 0.05.
    EXPECT_NEAR(cg::AlignedDistance(a, Substituted(a, 45), aligner, cigar), 45.0 / 900, 1e-12);
    // A copy called 60 bases longer at its start (another start codon) and 30 shorter at its end differs by its
    // substitutions only: the ends are free.
    std::string const longer = test::RandomSequence(60, rng) + Substituted(a, 18).substr(0, 870);
    EXPECT_NEAR(cg::AlignedDistance(a, longer, aligner, cigar), 18.0 / 870, 1e-12);
    // One base deleted inside: one column of 900 differs.
    std::string deleted = a;
    deleted.erase(450, 1);
    EXPECT_NEAR(cg::AlignedDistance(a, deleted, aligner, cigar), 1.0 / 900, 1e-12);
    // Unrelated sequences: beyond kMaxDistance (given up), or without enough overlap.
    double const far = cg::AlignedDistance(a, test::RandomSequence(900, rng), aligner, cigar);
    EXPECT_TRUE(std::isnan(far) || far >= 0.3) << far;
    EXPECT_TRUE(std::isnan(cg::AlignedDistance(a, "", aligner, cigar)));
}

TEST(CongenerGaps, AGeneOfAGenusGivesEachCopyItsNearestAndMedianCongener) {
    std::mt19937 rng(5);
    std::string const base = test::RandomSequence(1000, rng);
    // Taxa 1, 2 and 4 of genus 10: 2 differs from 1 at 10 positions, 4 at 30 others; taxon 3 alone in genus 11.
    std::string const one = base, two = Substituted(base, 10, 3, 7), four = Substituted(base, 30, 200, 11);
    std::vector<uint32_t> const taxids{ 1, 2, 3, 4 };
    std::vector<std::string> const seqs{ one, two, test::RandomSequence(1000, rng), four };
    auto genus = [](uint32_t taxid) { return taxid == 3 ? 11u : 10u; };
    std::vector<std::pair<uint32_t, cg::Gap>> out;
    cg::Stats stats;
    cg::ScanGene(7, taxids, seqs, genus, cg::Settings{}, 2, out, stats);
    std::map<uint32_t, cg::Gap> gaps(out.begin(), out.end());
    ASSERT_EQ(gaps.size(), 3u);  // not taxon 3: no congener
    EXPECT_EQ(gaps.at(1), (cg::Gap{ cg::Scaled(0.010), cg::Scaled(0.020), 2, 2 }));  // its nearest: taxon 2
    EXPECT_EQ(gaps.at(2), (cg::Gap{ cg::Scaled(0.010), cg::Scaled(0.025), 2, 1 }));  // 2 to 4: 40 differences
    EXPECT_EQ(gaps.at(4), (cg::Gap{ cg::Scaled(0.030), cg::Scaled(0.035), 2, 1 }));
    EXPECT_EQ(stats.genes, 1u);
    EXPECT_EQ(stats.alignments, 3u);
    EXPECT_EQ(stats.copies, 3u);
    EXPECT_EQ(stats.sampled_genera, 0u);
}

TEST(CongenerGaps, AGenusAboveTheLimitIsSampledAndStillFindsTheNearest) {
    std::mt19937 rng(7);
    std::string const base = test::RandomSequence(1200, rng);
    std::vector<uint32_t> taxids;
    std::vector<std::string> seqs;
    // Taxon 1 and its twin 2 (6 differences); taxa 3..12 each 60-120 differences from taxon 1, at other places.
    taxids.push_back(1);
    seqs.push_back(base);
    taxids.push_back(2);
    seqs.push_back(Substituted(base, 6, 2, 13));
    for (uint32_t t = 3; t <= 12; t++) {
        taxids.push_back(t);
        seqs.push_back(Substituted(base, 60 + 6 * t, 1 + t, 9));
    }
    cg::Settings settings;
    settings.all = 4;
    settings.nearest = 1;
    settings.sample = 3;
    std::vector<std::pair<uint32_t, cg::Gap>> out;
    cg::Stats stats;
    cg::ScanGene(2, taxids, seqs, [](uint32_t) { return 10u; }, settings, 3, out, stats);
    std::map<uint32_t, cg::Gap> gaps(out.begin(), out.end());
    EXPECT_EQ(stats.sampled_genera, 1u);
    ASSERT_TRUE(gaps.contains(1) && gaps.contains(2));
    EXPECT_EQ(gaps.at(1).min, cg::Scaled(6.0 / 1200));
    EXPECT_EQ(gaps.at(2).min, cg::Scaled(6.0 / 1200));
    EXPECT_EQ(gaps.at(1).nearest, 2u);  // the twins are each other's nearest, found by their sketches
    EXPECT_EQ(gaps.at(2).nearest, 1u);
    EXPECT_EQ(gaps.at(1).congeners, 11);
    // Each copy aligned against at most 1 + 3 others, so fewer pairs than all 66.
    EXPECT_LT(stats.alignments, 66u);
    // The same table for any number of threads.
    std::vector<std::pair<uint32_t, cg::Gap>> again;
    cg::Stats more;
    cg::ScanGene(2, taxids, seqs, [](uint32_t) { return 10u; }, settings, 1, again, more);
    EXPECT_EQ(again, out);
}

TEST(CongenerGaps, TheTableWritesAndReadsBack) {
    auto const table = cg::Table::FromRows({ { 5, 2, cg::Gap{ 120, 800, 3, 8 } }, { 2, 9, cg::Gap{ 0, 40, 1, 3 } },
                                             { 5, 1, cg::Gap{ 300, 300, 1, 6 } } });
    EXPECT_EQ(table.Copies(), 3u);
    EXPECT_EQ(table.Species(), 2u);
    std::ostringstream os;
    table.Write(os);
    EXPECT_EQ(os.str(), "taxid\tgene:min:median:congeners:nearest\n2\t9:0:40:1:3\n5\t1:300:300:1:6,2:120:800:3:8\n");
    cg::Table back;
    std::istringstream is(os.str());
    EXPECT_EQ(back.Read(is), "");
    ASSERT_TRUE(back.Find(5, 2).has_value());
    EXPECT_EQ(*back.Find(5, 2), (cg::Gap{ 120, 800, 3, 8 }));
    EXPECT_NEAR(back.Find(5, 2)->Median(), 0.08, 1e-12);
    EXPECT_FALSE(back.Find(5, 3).has_value());
    EXPECT_FALSE(back.Find(4, 1).has_value());
    EXPECT_FALSE(back.Find(99, 1).has_value());
    // A table of before 2026-10-08, without the nearest congener, reads with nearest 0; a trailing colon does not.
    cg::Table old;
    std::istringstream four("taxid\tgene:min:median:congeners\n7\t1:20:30:2\n");
    EXPECT_EQ(old.Read(four), "");
    EXPECT_EQ(*old.Find(7, 1), (cg::Gap{ 20, 30, 2, 0 }));
    std::istringstream trailing("7\t1:20:30:2:\n");
    EXPECT_NE(cg::Table().Read(trailing), "");
    std::istringstream bad("taxid\tgene:min:median:congeners\n7\t1:20:x:1\n");
    EXPECT_NE(cg::Table().Read(bad), "");
    std::istringstream range("7\t1:20000:30:1\n");
    EXPECT_NE(cg::Table().Read(range), "");
}

TEST(ForeignRates, TheTableWritesReadsBackAndSaturates) {
    auto const table = fr::Table::FromRows({ { 3, 4, fr::Rate{ 9, 4, 1 } }, { 3, 2, fr::Rate{ 2, 0, 0 } } });
    std::ostringstream os;
    table.Write(os);
    EXPECT_EQ(os.str(), "taxid\tgene:reads:foreign:foreign_genus\n3\t2:2:0:0,4:9:4:1\n");
    fr::Table back;
    std::istringstream is("# a comment\n" + os.str() + "8\t1:70000:69000:5\n");
    EXPECT_EQ(back.Read(is), "");
    ASSERT_NE(back.Find(3, 4), nullptr);
    EXPECT_NEAR(back.Find(3, 4)->ForeignShare(), 4.0 / 10, 1e-12);
    EXPECT_NEAR(back.Find(3, 4)->ForeignGenusShare(), 1.0 / 10, 1e-12);
    EXPECT_EQ(back.Find(8, 1)->reads, 65535);
    EXPECT_EQ(back.Find(8, 1)->foreign, 65535);
    std::istringstream inconsistent("3\t1:2:3:0\n");  // more foreign reads than reads
    EXPECT_NE(fr::Table().Read(inconsistent), "");
}

TEST(UntriedCandidates, TheyAreListedOnceInOrderApartFromTheFailedOnes) {
    std::vector<FailedCandidate> failed{ { 4, 2 } };
    std::vector<FailedCandidate> const attempted{ { 4, 2 }, { 7, 1 } };
    AddUntriedCandidates(failed, { 9, 7, 3, 9, 12 }, attempted, { 12 });
    ASSERT_EQ(failed.size(), 3u);
    EXPECT_EQ(FailedTag(failed), "4:2");
    EXPECT_EQ(UntriedTag(failed), "9,3");
    std::vector<FailedCandidate> many;
    std::vector<uint32_t> crowd;
    for (uint32_t t = 100; t < 120; t++) crowd.push_back(t);
    AddUntriedCandidates(many, crowd, {}, {});
    EXPECT_EQ(many.size(), kUntriedListed);
    // A read with untried candidates only gets no unmapped record (no ZF).
    SamEntry sam;
    EXPECT_FALSE(UnmappedRecord(sam, "r", many));
}

// A read of taxon 4's strain, ranked fourth by its seeds behind three congeners of the same genus: with --align_top 3 it
// aligns to a congener and its own species stays untried (ZC); the adaptive candidates try it and find it.
TEST(AdaptiveCandidates, ADivergentReadTriesItsCrowdsCongeners) {
    std::mt19937 rng(13);
    std::string const base = test::RandomSequence(600, rng);
    test::LoadedReference ref({ { 1, { Substituted(base, 30, 2, 19) } }, { 2, { Substituted(base, 30, 4, 19) } },
                                { 3, { Substituted(base, 30, 6, 19) } }, { 4, { base } }, { 5, { test::RandomSequence(600, rng) } } },
                              "adaptive");
    std::string const read = Substituted(base.substr(100, 150), 2, 40, 50);  // 2 differences from taxon 4
    std::string const rev = KmerUtils::ReverseComplement(read);
    auto anchors_of = [&] {
        AlignmentAnchorList anchors;
        for (uint32_t t : { 1u, 2u, 3u, 4u }) {
            ChainAlignmentAnchor anchor(t, 1, true);
            anchor.chain.emplace_back(100u, static_cast<uint16_t>(0), static_cast<uint16_t>(20));
            anchor.total_length = static_cast<uint16_t>(t < 4 ? 100 : 90);  // taxon 4 last, in the crowd (>= 0.8)
            anchors.push_back(anchor);
        }
        ChainAlignmentAnchor other(5, 1, true);  // another genus, far behind
        other.chain.emplace_back(100u, static_cast<uint16_t>(0), static_cast<uint16_t>(20));
        other.total_length = 30;
        anchors.push_back(other);
        return anchors;
    };
    WFA2Wrapper2 aligner(4, 6, 2, 1000);
    auto genera = std::make_shared<std::vector<uint32_t> const>(std::vector<uint32_t>{ 0, 10, 10, 10, 10, 11 });
    std::string id = "r";
    for (size_t extra : { 0u, 7u }) {
        SCOPED_TRACE(extra);
        SimpleAlignmentHandler handler(*ref.loader, aligner, 31, 3, 0.9, false);
        handler.SetAnchoredAlignment(false);
        handler.SetAdaptiveCandidates(extra, genera);
        auto anchors = anchors_of();
        AlignmentResultList results;
        handler(anchors, results, read, rev, 3, id);
        ASSERT_FALSE(results.empty());
        EXPECT_EQ(handler.Crowding(), 4);
        if (extra == 0) {
            EXPECT_NE(results.front().Taxid(), 4u);
            EXPECT_EQ(handler.Untried(), (std::vector<uint32_t>{ 4 }));
            EXPECT_EQ(handler.m_adaptive_alignments, 0u);
        } else {
            EXPECT_EQ(results.front().Taxid(), 4u);
            EXPECT_TRUE(handler.Untried().empty());
            EXPECT_EQ(handler.m_adaptive_alignments, 1u);  // taxon 5 is of another genus and not in the crowd
        }
    }
    // A read of taxon 1 itself (no divergence) tries nothing more.
    SimpleAlignmentHandler handler(*ref.loader, aligner, 31, 3, 0.9, false);
    handler.SetAnchoredAlignment(false);
    handler.SetAdaptiveCandidates(7, genera);
    std::string const own = ref.Gene(1, 1).substr(100, 150);
    auto anchors = anchors_of();
    AlignmentResultList results;
    handler(anchors, results, own, KmerUtils::ReverseComplement(own), 3, id);
    EXPECT_EQ(results.front().Taxid(), 1u);
    EXPECT_EQ(handler.m_adaptive_alignments, 0u);
}

TEST(CopyFeatures, TheGapsForeignRatesAndUntriedCandidatesOfATaxon) {
    std::mt19937 rng(17);
    test::LoadedReference ref({ { 1, { test::RandomSequence(600, rng), test::RandomSequence(600, rng) } },
                                { 2, { test::RandomSequence(600, rng) } }, { 3, { test::RandomSequence(600, rng) } } },
                              "copy features");
    std::string sam = ref.Header();
    // Taxon 1: an exact read on gene 1, one with 12 mismatches of 150 (0.08), and one on gene 2, which has no gap.
    sam += Record(ref, "a", 1, 1, 0, 150);
    sam += Record(ref, "b", 1, 1, 200, 150, "69=12X69=", "ZA:Z:*\tZC:Z:2,3");
    sam += Record(ref, "c", 1, 2, 0, 150);
    sam += Record(ref, "d", 2, 1, 0, 150);
    auto const without = Features(Profile(ref, sam).GetTaxa().at(1));
    for (auto const* name : { "gap_informative_share", "gap_within_min_share", "gap_within_median_share", "gap_position",
                              "foreign_scanned_share", "foreign_copy_share", "foreign_genus_copy_share" }) {
        EXPECT_EQ(without.at(name), -1.0) << name;
    }
    // Gene 1 of taxon 1: nearest congener 0.05 away, median 0.10; scanned with 10 reads, 5 of other species, 1 of another genus.
    ref.loader->SetCongenerGaps(cg::Table::FromRows({ { 1, 1, cg::Gap{ 500, 1000, 4 } }, { 2, 1, cg::Gap{ 20, 30, 1 } } }));
    ref.loader->SetForeignRates(fr::Table::FromRows({ { 1, 1, fr::Rate{ 10, 5, 1 } } }));
    for (size_t threads : { 1, 3 }) {
        SCOPED_TRACE(threads);
        auto const profile = Profile(ref, sam, threads);
        auto const f1 = Features(profile.GetTaxa().at(1));
        auto const f2 = Features(profile.GetTaxa().at(2));
        EXPECT_NEAR(f1.at("gap_informative_share"), 2.0 / 3, 1e-12);
        EXPECT_NEAR(f1.at("gap_within_min_share"), 0.5, 1e-12);     // 0 < 0.05, 0.08 not
        EXPECT_NEAR(f1.at("gap_within_median_share"), 1.0, 1e-12);  // both below 0.10
        EXPECT_NEAR(f1.at("gap_position"), 0.05, 1e-12);            // the lower median: 0 / 0.05, bin 0
        EXPECT_NEAR(f1.at("foreign_scanned_share"), 2.0 / 3, 1e-12);
        EXPECT_NEAR(f1.at("foreign_copy_share"), 5.0 / 11, 1e-9);
        EXPECT_NEAR(f1.at("foreign_genus_copy_share"), 1.0 / 11, 1e-9);
        // Taxon 2's copy is too near its congener (0.002) to tell: no informative record; read b never tried it (ZC).
        EXPECT_EQ(f2.at("gap_informative_share"), 0.0);
        EXPECT_EQ(f2.at("gap_within_min_share"), -1.0);
        EXPECT_EQ(f2.at("foreign_scanned_share"), 0.0);
        EXPECT_EQ(f2.at("foreign_copy_share"), -1.0);
        EXPECT_NEAR(f2.at("untried_candidate_rate"), 0.5, 1e-12);
        EXPECT_EQ(f1.at("untried_candidate_rate"), 0.0);
        EXPECT_FALSE(profile.GetTaxa().contains(3));  // untried alone gives no taxon
    }
}
