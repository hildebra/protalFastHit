// Unit tests for the strain alleles (2026-10-08, docs/claude/2026-10-08-strain-alleles): the genomes that may give them
// (the accession hash, the same in Python), the table and its file, an allele's edits from its alignment against the
// representative, the build's sample and choice, what an allele explains of a read, the alignment score's shift, the
// handler's scores, and the profiler's alleles features. And the polymorphic sites (2026-10-09,
// docs/claude/2026-10-09-polymorphic-sites): a copy's sites, a read at them, the site shift, the unsure reads the
// handler settles with it, and the profiler's polymorphic features.
#include <gtest/gtest.h>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include "Core/AlignmentStrategy.h"
#include "Profiling/Profiler.h"
#include "SequenceUtils/KmerIterator.h"
#include "SequenceUtils/Minimizer.h"
#include "SequenceUtils/StrainAllelesBuild.h"
#include "TestReference.h"

using namespace protal;
namespace sa = protal::strain_alleles;

namespace {
    // seq with `n` substitutions at evenly spaced positions from `first` (every `step`-th base).
    std::string Substituted(std::string seq, size_t n, size_t first, size_t step) {
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

    profiler::MicrobialProfile Profile(test::LoadedReference const& ref, std::string const& sam, size_t threads = 1) {
        profiler::Profiler profiler(*ref.loader);
        profiler.SetDepthIdentityMargin(0.08);
        profiler::MicrobialProfile profile(*ref.loader);
        auto const path = ref.Write("alleles" + std::to_string(threads) + ".sam", sam);
        std::ostringstream rejected;
        std::string const error = profiler.ProfileSam(path, profile, std::optional<std::reference_wrapper<std::ostream>>{ rejected },
                                                      2, 2, 0.0, 15, 0, false, threads);
        EXPECT_EQ(error, "");
        return profile;
    }

    std::string SamRecord(std::string const& qname, uint32_t taxid, uint32_t gene, int start, std::string const& cigar,
                          std::string const& seq) {
        return qname + "\t0\t" + std::to_string(taxid) + "_" + std::to_string(gene) + '\t' + std::to_string(start + 1) + "\t60\t" +
               cigar + "\t*\t0\t0\t" + seq + '\t' + std::string(seq.size(), 'I') + "\tZA:Z:*\n";
    }

    char Other(char base) { return base == 'A' ? 'C' : 'A'; }
}

// The genomes that may give alleles: a hash of the accession, the same in scripts/mini_db/gtdb_to_protal_db.py
// (allele_genome; scripts/mini_db/test_gtdb_build.py pins the same values).
TEST(StrainAlleles, TheGenomeSplitIsTheHashOfTheAccession) {
    EXPECT_EQ(sa::Fnv1a("GCF_000005845.2"), 5130540707165718432ull);
    EXPECT_EQ(sa::Fnv1a("GCA_900000000.1"), 12847648632108725385ull);
    EXPECT_EQ(sa::Fnv1a(""), 1469598103934665603ull);
    EXPECT_TRUE(sa::AlleleGenome("GCF_000005845.2", 0.5));   // 0.278
    EXPECT_TRUE(sa::AlleleGenome("GCA_000001405.29", 0.5));  // 0.090
    EXPECT_FALSE(sa::AlleleGenome("GCA_900000000.1", 0.5));  // 0.696
    EXPECT_FALSE(sa::AlleleGenome("GCF_000195955.2", 0.5));  // 0.508
    EXPECT_TRUE(sa::AlleleGenome("GCF_000195955.2", 0.51));
    EXPECT_TRUE(sa::AlleleGenome("GCA_900000000.1", 1));
    EXPECT_FALSE(sa::AlleleGenome("GCF_000005845.2", 0));
}

TEST(StrainAlleles, TheTableWritesAndReadsBack) {
    sa::Allele a{ 0, 900, { sa::Edit::Substitution(12, 2), sa::Edit::Indel(sa::kInsertion, 300, 3), sa::Edit::Indel(sa::kDeletion, 450, 6) } };
    sa::Allele b{ 20, 880, { sa::Edit::Substitution(25, 3) } };
    sa::Allele c{ 0, 1000, { sa::Edit::Substitution(999, 0) } };
    auto const table = sa::Table::FromRows({ { 5, 2, { a, b } }, { 3, 7, { c } }, { 5, 1, { b } } });
    EXPECT_EQ(table.Copies(), 3u);
    EXPECT_EQ(table.Species(), 2u);
    EXPECT_EQ(table.Alleles(), 4u);
    EXPECT_EQ(table.Edits(), 6u);
    ASSERT_EQ(table.Of(5, 2).size(), 2u);
    EXPECT_EQ(table.Get(table.Of(5, 2)[0]).edits[1], sa::Edit::Indel(sa::kInsertion, 300, 3));
    EXPECT_TRUE(table.Of(5, 3).empty());
    EXPECT_TRUE(table.Of(9, 1).empty());
    std::ostringstream os;
    table.Write(os);
    std::string const text = os.str();
    EXPECT_NE(text.find("\n3\t7\t0-1000:999A\n5\t1\t20-880:25T\n5\t2\t0-900:12G,300i3,450d6;20-880:25T\n"), std::string::npos) << text;
    sa::Table back;
    std::istringstream is(text);
    ASSERT_EQ(back.Read(is), "");
    std::ostringstream again;
    back.Write(again);
    EXPECT_EQ(again.str(), text);
    for (std::string const bad : { "5\t2\t0-900:12G,", "5\t2\t900-100:12G", "5\t2\t0-900:30G,12A", "5\t2\t0-900:12N", "5\t2\t",
                                   "5\t2\t0-900:12x3", "5\t0-900:1A" }) {
        sa::Table t;
        std::istringstream bad_is(bad + std::string("\n"));
        EXPECT_NE(t.Read(bad_is), "") << bad;
    }
}

// A copy's allele against the representative's copy: a substitution, an insertion and a deletion where the copy covers
// the representative, its overhang and the representative's uncovered ends left out.
TEST(StrainAlleles, AnAllelesEditsComeFromItsAlignmentToTheRepresentative) {
    std::mt19937 rng(23);
    std::string const rep = test::RandomSequence(800, rng);
    std::string copy = rep.substr(30, 740);  // covers [30, 770)
    copy[100 - 30] = Other(copy[100 - 30]);
    copy.erase(400 - 30, 2);                 // the representative's 400, 401 deleted
    copy.insert(300 - 30, "GTC");            // three bases before the representative's 300
    copy = test::RandomSequence(12, rng) + copy + test::RandomSequence(7, rng);  // an overhang at both ends
    WFA2Wrapper2 aligner(4, 6, 2, 0);
    std::string ops;
    double divergence = 0;
    auto const allele = sa::Diff(rep, copy, aligner, ops, divergence);
    ASSERT_TRUE(allele);
    EXPECT_GE(allele->begin, 28u);  // about where the copy starts: free ends may take a base or two
    EXPECT_LE(allele->begin, 31u);
    EXPECT_GE(allele->end, 769u);
    EXPECT_LE(allele->end, 772u);
    // The indels where WFA2 puts them: in a repeat, a base or two apart.
    ASSERT_EQ(allele->edits.size(), 3u);
    EXPECT_EQ(allele->edits[0], sa::Edit::Substitution(100, sa::BaseCode(copy[12 + 70])));
    EXPECT_EQ(allele->edits[1].GetKind(), sa::kInsertion);
    EXPECT_EQ(allele->edits[1].Length(), 3u);
    EXPECT_NEAR(allele->edits[1].pos, 300, 3);
    EXPECT_EQ(allele->edits[2].GetKind(), sa::kDeletion);
    EXPECT_EQ(allele->edits[2].Length(), 2u);
    EXPECT_NEAR(allele->edits[2].pos, 400, 3);
    EXPECT_NEAR(divergence, 6.0 / 743, 0.0005);
    // An unrelated copy is no allele; nor is one that covers too little.
    EXPECT_FALSE(sa::Diff(rep, test::RandomSequence(800, rng), aligner, ops, divergence));
    EXPECT_FALSE(sa::Diff(rep, rep.substr(0, 300), aligner, ops, divergence));
    // An insertion right before a substituted base: both at one position, the substitution first as the table keeps
    // them (the alignment gives the insertion first), so that the table reads back (an r226-like pipeline run had one).
    std::string both = rep;
    both[500] = Other(both[500]);
    both.insert(500, "TT");
    auto const pair = sa::Diff(rep, both, aligner, ops, divergence);
    ASSERT_TRUE(pair);
    EXPECT_TRUE(std::is_sorted(pair->edits.begin(), pair->edits.end()));
    std::ostringstream os;
    sa::Table::FromRows({ { 1, 1, { *pair } } }).Write(os);
    sa::Table back;
    std::istringstream is(os.str());
    EXPECT_EQ(back.Read(is), "") << os.str();
    EXPECT_EQ(back.Alleles(), 1u);
}

// The index's seeds of a copy's deep alleles (AlleleSeeds, --index_alleles): k-mers of the representative's copy with
// the allele's substitutions that the copy does not yield, none over an indel, none the representative or an earlier
// allele gave, and none of an allele below the divergence asked for.
TEST(StrainAlleles, DeepAllelesGiveTheIndexSeedsTheRepresentativeLacks) {
    using E = sa::Edit;
    constexpr size_t kK = 31, kM = 15;
    ClosedSyncmer syncmer{ kM, 7, 2, true };
    SimpleKmerHandler<ClosedSyncmer> handler{ kK, kM, syncmer };
    std::mt19937 rng(21);
    std::string const rep = test::RandomSequence(600, rng);
    KmerList own;
    handler(std::string_view(rep), own);
    ASSERT_GT(own.size(), 20u);
    std::vector<uint64_t> own_keys;
    for (auto const& p : own) own_keys.push_back(p.first);
    std::sort(own_keys.begin(), own_keys.end());
    auto other = [&rep](size_t p) { return static_cast<uint8_t>((sa::BaseCode(rep[p]) + 1) & 3); };
    // deep: 12 substitutions in 600 (2%); near: 2 (0.3%); indel: 12 substitutions and an insertion before 300.
    sa::Allele deep{ 0, 600, {} }, near{ 0, 600, { E::Substitution(123, other(123)), E::Substitution(437, other(437)) } }, indel{ 0, 600, {} };
    for (uint16_t p = 50; p < 600; p += 50) {
        deep.edits.push_back(E::Substitution(p, other(p)));
        indel.edits.push_back(E::Substitution(p, other(p)));
    }
    indel.edits.push_back(E::Indel(sa::kInsertion, 300, 3));
    std::sort(indel.edits.begin(), indel.edits.end());
    auto const table = sa::Table::FromRows({ { 1, 1, { deep, near } }, { 1, 2, { indel } }, { 2, 1, { near } } });
    sa::SeedScratch scratch;
    KmerList seeds;
    // The deep allele's seeds: new k-mers, each window holding one of its substitutions, at the representative's positions.
    EXPECT_EQ(sa::AlleleSeeds(table, 0.01, 1, 1, rep, own, kK, handler, scratch, seeds), 1u);
    ASSERT_GT(seeds.size(), 10u);
    std::set<uint64_t> distinct;
    for (auto const& [key, pos] : seeds) {
        EXPECT_FALSE(std::binary_search(own_keys.begin(), own_keys.end(), static_cast<uint64_t>(key)));
        EXPECT_TRUE(distinct.insert(key).second);
        EXPECT_LE(pos + kK, 600u);
        bool holds = false;
        for (auto const& e : deep.edits) holds = holds || (e.pos >= pos && e.pos < pos + kK);
        EXPECT_TRUE(holds) << pos;
    }
    // Below 0.01 the near allele adds nothing; at 0 it does, and only k-mers the deep one did not give.
    seeds.clear();
    EXPECT_EQ(sa::AlleleSeeds(table, 0.0, 1, 1, rep, own, kK, handler, scratch, seeds), 2u);
    EXPECT_GT(seeds.size(), distinct.size());
    std::set<uint64_t> all;
    for (auto const& [key, pos] : seeds) EXPECT_TRUE(all.insert(key).second);
    // No window over the insertion (bases 299 and 300 are apart in the allele).
    seeds.clear();
    EXPECT_EQ(sa::AlleleSeeds(table, 0.01, 1, 2, rep, own, kK, handler, scratch, seeds), 1u);
    ASSERT_GT(seeds.size(), 5u);
    for (auto const& [key, pos] : seeds) EXPECT_FALSE(pos < 300 && 300 < pos + kK) << pos;
    // A copy without alleles, or without a deep one: nothing.
    seeds.clear();
    EXPECT_EQ(sa::AlleleSeeds(table, 0.01, 2, 1, rep, own, kK, handler, scratch, seeds), 0u);
    EXPECT_EQ(sa::AlleleSeeds(table, 0.01, 3, 1, rep, own, kK, handler, scratch, seeds), 0u);
    EXPECT_TRUE(seeds.empty());
    EXPECT_NEAR(sa::Divergence(table.Get(table.Of(1, 1)[0])), 11.0 / 600, 1e-9);
}

// The build's sample of a copy's alleles (the least hashes) and its choice (the sample's coverage, greedily) do not
// depend on the order the alleles come in; an allele as far as the nearest congener's copy, of a genome outside the
// share, identical to the representative or seen before is not offered.
TEST(StrainAlleles, TheBuildsSampleAndChoiceDoNotDependOnTheOrder) {
    using E = sa::Edit;
    sa::Allele const a{ 0, 100, { E::Substitution(10, 0), E::Substitution(20, 1), E::Substitution(30, 2) } };  // 3 from the rep
    sa::Allele const b{ 0, 100, { E::Substitution(10, 0) } };                                                   // 1, near a
    sa::Allele const c{ 0, 100, { E::Substitution(40, 3), E::Substitution(50, 3) } };                           // 2, far from a
    // Stored, b saves a read of a one of its three edits (a shares b's one edit, more than half of b's), a saves a read
    // of b nothing (b lacks two of a's three: the read keeps the representative), c nothing of either.
    EXPECT_EQ(sa::Coverage(a, b), sa::kCoverageUnit / 3);
    EXPECT_EQ(sa::Coverage(b, a), 0u);
    EXPECT_EQ(sa::Coverage(a, a), sa::kCoverageUnit);
    EXPECT_EQ(sa::Coverage(c, a), 0u);
    // So b covers 1 + 1/3 of the sample, a and c 1 each: b first, then c (1 against the 2/3 a still adds), then a.
    EXPECT_EQ(sa::Select({ { 1, a }, { 2, b }, { 3, c } }, 2), (std::vector<sa::Allele>{ b, c }));
    EXPECT_EQ(sa::Select({ { 1, a }, { 2, b }, { 3, c } }, 5), (std::vector<sa::Allele>{ b, c, a }));
    // Farthest first (until 2026-10-09) took a and c: a deep lineage's alleles over a strain's near relatives.
    sa::Allele const a2 = a;
    EXPECT_EQ(sa::Select({ { 1, a }, { 2, a2 }, { 3, c } }, 5), (std::vector<sa::Allele>{ a, c }));  // a twin adds nothing
    std::vector<sa::Allele> many;
    for (uint16_t p = 1; p <= 12; p++) many.push_back(sa::Allele{ 0, 200, { E::Substitution(p, 0), E::Substitution(static_cast<uint16_t>(p + 100), 1) } });
    auto rows_of = [&](std::vector<sa::Allele> order, size_t threads) {
        sa::Collector collector(4);
        std::vector<std::thread> workers;
        for (size_t t = 0; t < threads; t++) {
            workers.emplace_back([&, t] {
                for (size_t i = t; i < order.size(); i += threads) collector.Offer(sa::Collector::Key(7, 3), order[i]);
            });
        }
        for (auto& w : workers) w.join();
        return collector.Rows(2);
    };
    auto const forward = rows_of(many, 1);
    std::reverse(many.begin(), many.end());
    EXPECT_EQ(rows_of(many, 1), forward);
    EXPECT_EQ(rows_of(many, 3), forward);
    ASSERT_EQ(forward.size(), 1u);
    EXPECT_EQ(std::get<2>(forward.front()).size(), 2u);
    // AddRecord's filters.
    std::mt19937 rng(5);
    std::string const rep = test::RandomSequence(500, rng);
    std::string const strain = Substituted(rep, 5, 50, 80);  // 0.01 from rep
    WFA2Wrapper2 aligner(4, 6, 2, 0);
    std::string ops;
    sa::Collector collector(4);
    sa::Stats stats;
    sa::Settings settings;
    settings.share = 0.5;
    sa::AddRecord(collector, stats, settings, 1, 1, "GCA_900000000.1", rep, strain, -1, aligner, ops);  // outside the share
    sa::AddRecord(collector, stats, settings, 1, 1, "GCF_000005845.2", rep, rep, -1, aligner, ops);     // identical
    sa::AddRecord(collector, stats, settings, 1, 1, "GCF_000005845.2", rep, strain, 0.005, aligner, ops);  // beyond the congener
    sa::AddRecord(collector, stats, settings, 1, 2, "GCF_000005845.2", rep, strain, 0.05, aligner, ops);   // kept
    sa::AddRecord(collector, stats, settings, 1, 2, "GCA_000001405.29", rep, strain, 0.05, aligner, ops);  // the same sequence again
    EXPECT_EQ(stats.records, 5u);
    EXPECT_EQ(stats.outside_share, 1u);
    EXPECT_EQ(stats.identical, 1u);
    EXPECT_EQ(stats.beyond_congener, 1u);
    EXPECT_EQ(stats.offered, 1u);
    EXPECT_EQ(stats.duplicates, 1u);
    auto const rows = collector.Rows(4);
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_EQ(std::get<1>(rows.front()), 2u);
    EXPECT_EQ(std::get<2>(rows.front()).front().edits.size(), 5u);
}

// What an allele makes of a read's differences: those it explains (its substitution with the read's base, its indel of
// the read's kind and length nearby) and those of its own the read lacks; a read's differences from a SAM record and from
// an alignment's columns are the same.
TEST(StrainAlleles, AnAlleleExplainsAReadsDifferences) {
    using E = sa::Edit;
    auto const table = sa::Table::FromRows({ { 4, 1, { sa::Allele{ 0, 600, { E::Substitution(110, 2), E::Indel(sa::kDeletion, 150, 2),
                                                                               E::Substitution(400, 1) } } } } });
    // A read on [100, 250): X at 110 (G, the allele's), a deletion of 2 at 152 (the allele's at 150, placed apart), X at 201.
    std::string seq(150 - 2, 'A');
    seq[10] = 'G';
    seq[99] = 'T';
    sa::ReadDiffs diffs;
    ASSERT_TRUE(sa::FromSamRecord("10=1X41=2D47=1X48=", 101, seq, diffs));
    EXPECT_EQ(diffs.begin, 100u);
    EXPECT_EQ(diffs.end, 250u);
    EXPECT_EQ(diffs.columns, 150u);
    ASSERT_EQ(diffs.diffs.size(), 3u);
    auto const best = sa::BestAllele(table, 4, 1, diffs);
    EXPECT_EQ(best.allele, 0);
    EXPECT_EQ(best.explained, 2u);
    EXPECT_EQ(best.shift, -2);  // 400 lies beyond the read
    // The same from the alignment's columns, the read as aligned.
    sa::ReadDiffs columns;
    std::string ops = std::string(10, 'M') + "X" + std::string(41, 'M') + "DD" + std::string(47, 'M') + "X" + std::string(48, 'M');
    sa::FromColumns(ops, seq, 100, columns);
    EXPECT_EQ(columns.diffs, diffs.diffs);
    EXPECT_EQ(columns.columns, diffs.columns);
    // A read with the reference's base at 110 and no deletion: the allele would add two differences; the reference stays.
    std::string const plain(150, 'A');
    ASSERT_TRUE(sa::FromSamRecord("150=", 101, plain, diffs));
    auto const none = sa::BestAllele(table, 4, 1, diffs);
    EXPECT_EQ(none.allele, -1);
    EXPECT_EQ(none.shift, 0);
    EXPECT_EQ(sa::Explain(table.Get(table.Of(4, 1)[0]), diffs).contradicted, 2u);
    // Another base at 110 is no match; a copy without alleles gives nothing.
    seq[10] = 'C';
    ASSERT_TRUE(sa::FromSamRecord("10=1X139=", 101, seq.substr(0, 150 - 2) + "AA", diffs));
    EXPECT_EQ(sa::BestAllele(table, 4, 1, diffs).explained, 0u);
    EXPECT_EQ(sa::BestAllele(table, 4, 2, diffs).allele, -1);
}

// The site shift counts half differences: two explained differences (-4) are two matches, one at a polymorphic site (-1)
// half of one.
TEST(StrainAlleles, TheScoreCountsTheShiftAsMatches) {
    AlignmentInfo info;
    info.matches = 97;
    info.mismatches = 3;
    EXPECT_EQ(info.Score(), -12);
    info.site_shift = -4;
    EXPECT_EQ(info.Score(), -4);
    EXPECT_EQ(info.Score(2, 3, 1, 2), 97 * 2 - 3 * 3 + 2 * 5);
    info.site_shift = -1;
    EXPECT_EQ(info.Score(), -10);
    info.site_shift_pending = -3;
    info.Reset();
    EXPECT_EQ(info.site_shift, 0);
    EXPECT_EQ(info.site_shift_pending, 0);
}

// A copy's polymorphic sites, how a read stands at them, and a candidate's site shift. Two alleles of copy (4, 1): A1 over
// the whole copy with G at 110 and T at 200; A2 over [100, 400) with C at 110, a deletion of 2 at 300 and A at 350.
TEST(StrainAlleles, ThePolymorphicSitesOfACopy) {
    using E = sa::Edit;
    auto const table = sa::Table::FromRows({ { 4, 1, { sa::Allele{ 0, 600, { E::Substitution(110, sa::BaseCode('G')),
                                                                              E::Substitution(200, sa::BaseCode('T')) } },
                                                         sa::Allele{ 100, 400, { E::Substitution(110, sa::BaseCode('C')),
                                                                                 E::Indel(sa::kDeletion, 300, 2),
                                                                                 E::Substitution(350, sa::BaseCode('A')) } } } } });
    auto const alleles = table.Of(4, 1);
    sa::Polymorphism poly;
    poly.Set(table, alleles, 100, 250);
    ASSERT_EQ(poly.Sites().size(), 2u);
    EXPECT_EQ(poly.Sites()[0].pos, 110u);
    EXPECT_EQ(poly.Sites()[0].bases, (1u << sa::BaseCode('G')) | (1u << sa::BaseCode('C')));
    EXPECT_EQ(poly.Sites()[1].pos, 200u);
    EXPECT_EQ(poly.Find(150), nullptr);
    EXPECT_EQ(poly.Cover(150), 2u);
    EXPECT_EQ(poly.Cover(50), 1u);
    EXPECT_EQ(poly.Cover(650), 0u);
    poly.Set(table, alleles, 280, 330);
    ASSERT_EQ(poly.Sites().size(), 2u);  // the deletion's two bases
    EXPECT_TRUE(poly.Sites()[0].indel && poly.Sites()[1].indel);
    EXPECT_TRUE(poly.IndelNear(304));
    EXPECT_FALSE(poly.IndelNear(310));

    // A read on [100, 250): C at 110 (A2's), A at 200 (no allele's: the species varies there), G at 150 (fixed).
    std::string seq(150, 'A');
    seq[10] = 'C';
    seq[50] = 'G';
    seq[100] = 'A';
    sa::ReadDiffs diffs;
    ASSERT_TRUE(sa::FromSamRecord("10=1X39=1X49=1X49=", 101, seq, diffs));
    poly.Set(table, alleles, diffs.begin, diffs.end);
    auto c = sa::CountSites(poly, diffs);
    EXPECT_EQ(c.sites, 2u);
    EXPECT_EQ(c.known, 1u);
    EXPECT_EQ(c.novel, 1u);
    // Its shift: A2 explains 110 (-1; A1 would contradict both), 200 is variable (half), 150 fixed: -3 half differences.
    std::vector<uint8_t> used;
    auto s = sa::ShiftOf(table, 4, 1, diffs, poly, used);
    EXPECT_EQ(s.best, -1);
    EXPECT_EQ(s.known, 0u);
    EXPECT_EQ(s.variable, 1u);
    EXPECT_EQ(s.Half(), -3);
    // C at 110 and T at 200: A2 explains 110, A1 has the T (a known variant though A1 contradicts 110): two matches.
    seq[100] = 'T';
    ASSERT_TRUE(sa::FromSamRecord("10=1X89=1X49=", 101, seq.substr(0, 50) + std::string(1, 'A') + seq.substr(51), diffs));
    s = sa::ShiftOf(table, 4, 1, diffs, poly, used);
    EXPECT_EQ(s.best, -1);
    EXPECT_EQ(s.known, 1u);
    EXPECT_EQ(s.Half(), -4);
    poly.Set(table, alleles, diffs.begin, diffs.end);
    c = sa::CountSites(poly, diffs);
    EXPECT_EQ(c.sites, 2u);
    EXPECT_EQ(c.known, 2u);
    // The representative's bases: both sites covered, neither way; no shift.
    ASSERT_TRUE(sa::FromSamRecord("150=", 101, std::string(150, 'A'), diffs));
    c = sa::CountSites(poly, diffs);
    EXPECT_EQ(c.sites, 2u);
    EXPECT_EQ(c.known + c.novel, 0u);
    EXPECT_EQ(sa::ShiftOf(table, 4, 1, diffs, poly, used).Half(), 0);
    // A read on [280, 330) with A2's deletion (placed at 301): known at both deleted bases; a deletion of 1 at 300: the
    // site at 300 is under it (an indel where an allele has one: known), and 301 within reach of it too.
    ASSERT_TRUE(sa::FromSamRecord("21=2D27=", 281, std::string(48, 'A'), diffs));
    poly.Set(table, alleles, diffs.begin, diffs.end);
    c = sa::CountSites(poly, diffs);
    EXPECT_EQ(c.sites, 2u);
    EXPECT_EQ(c.known, 2u);
    s = sa::ShiftOf(table, 4, 1, diffs, poly, used);
    EXPECT_EQ(s.Half(), -2);  // A2's deletion, within its tolerance: explained
    // A substitution where the species is fixed and no allele covers: nothing (copy 4, 2 has no alleles).
    EXPECT_EQ(sa::ShiftOf(table, 4, 2, diffs, poly, used).Half(), 0);
}

// A read of a known strain of taxon 1 (its allele: the gene before taxon 1's own six differences) fits taxon 2 better by
// the references alone (four differences against six), within kUnsureMismatches: unsure, so the candidates take their
// site shifts, and it scores on taxon 1 as on the strain's gene.
TEST(StrainAlleles, AReadOfAKnownStrainScoresOnItsSpecies) {
    std::mt19937 rng(31);
    std::string const base = test::RandomSequence(600, rng);
    std::string const rep1 = Substituted(base, 6, 105, 20);
    std::string const rep2 = Substituted(base, 4, 110, 30);
    test::LoadedReference ref({ { 1, { rep1 } }, { 2, { rep2 } } }, "allele scores");
    std::vector<sa::Edit> edits;
    for (size_t k = 0; k < 6; k++) {
        size_t const p = 105 + k * 20;
        edits.push_back(sa::Edit::Substitution(static_cast<uint16_t>(p), sa::BaseCode(base[p])));
    }
    ref.loader->SetStrainAlleles(sa::Table::FromRows({ { 1, 1, { sa::Allele{ 0, 600, edits } } } }));
    std::string const read = base.substr(100, 150);
    std::string const rev = KmerUtils::ReverseComplement(read);
    WFA2Wrapper2 aligner(4, 6, 2, 1000);
    std::string id = "r";
    for (bool alleles : { false, true }) {
        SCOPED_TRACE(alleles);
        SimpleAlignmentHandler handler(*ref.loader, aligner, 31, 3, 0.9, false);
        handler.SetAnchoredAlignment(false);
        handler.SetAlleleScores(alleles);
        AlignmentAnchorList anchors;
        for (uint32_t t : { 2u, 1u }) {
            ChainAlignmentAnchor anchor(t, 1, true);
            anchor.chain.emplace_back(100u, static_cast<uint16_t>(0), static_cast<uint16_t>(20));
            anchor.total_length = 100;
            anchors.push_back(anchor);
        }
        AlignmentResultList results;
        handler(anchors, results, read, rev, 3, id);
        ASSERT_EQ(results.size(), 2u);
        EXPECT_EQ(results.front().Taxid(), alleles ? 1u : 2u);
        auto const& on1 = results.front().Taxid() == 1 ? results.front() : results.back();
        EXPECT_EQ(on1.GetAlignmentInfo().mismatches, 6);  // its CIGAR stays the reference's
        EXPECT_EQ(on1.GetAlignmentInfo().site_shift, alleles ? -12 : 0);
        EXPECT_EQ(handler.m_allele_scored, alleles ? 1u : 0u);
        EXPECT_EQ(handler.m_allele_shifted, alleles ? 1u : 0u);
        EXPECT_EQ(handler.m_unsure_reads, alleles ? 1u : 0u);
        EXPECT_EQ(handler.m_settled_reads, alleles ? 1u : 0u);
    }
    // A read of taxon 2's own gene: ten differences on taxon 1, none on taxon 2. Taxon 1's allele would explain six of
    // them, but the read is clear (beyond kUnsureMismatches): no candidate takes its shift.
    std::string const own = rep2.substr(100, 150);
    SimpleAlignmentHandler handler(*ref.loader, aligner, 31, 3, 0.5, false);
    handler.SetAnchoredAlignment(false);
    handler.SetAlleleScores(true);
    AlignmentAnchorList anchors;
    for (uint32_t t : { 2u, 1u }) {
        ChainAlignmentAnchor anchor(t, 1, true);
        anchor.chain.emplace_back(100u, static_cast<uint16_t>(0), static_cast<uint16_t>(20));
        anchor.total_length = 100;
        anchors.push_back(anchor);
    }
    AlignmentResultList results;
    handler(anchors, results, own, KmerUtils::ReverseComplement(own), 3, id);
    ASSERT_EQ(results.size(), 2u);
    EXPECT_EQ(results.front().Taxid(), 2u);
    auto const& on1 = results.back();
    EXPECT_EQ(on1.GetAlignmentInfo().mismatches, 10);
    EXPECT_EQ(on1.GetAlignmentInfo().site_shift_pending, -12);
    EXPECT_EQ(on1.GetAlignmentInfo().site_shift, 0);
    EXPECT_EQ(handler.m_unsure_reads, 0u);
}

// The profiler's alleles features: the share of a taxon's kept records on copies with alleles (no feature group), the share
// of their differences the best allele explains, and the identity it gains; -1 without the table, 0 without alleles.
TEST(StrainAlleles, TheAllelesFeaturesOfATaxon) {
    std::mt19937 rng(41);
    test::LoadedReference ref({ { 1, { test::RandomSequence(600, rng) } }, { 2, { test::RandomSequence(600, rng) } } }, "allele features");
    std::string const gene = ref.Gene(1, 1);
    // Taxon 1's allele differs at 10 and 20 (by other bases).
    auto const table = sa::Table::FromRows({ { 1, 1, { sa::Allele{ 0, 600, { sa::Edit::Substitution(10, sa::BaseCode(Other(gene[10]))),
                                                                              sa::Edit::Substitution(20, sa::BaseCode(Other(gene[20]))) } } } } });
    std::string sam = ref.Header();
    // a: the allele's two bases (explained); b: a difference at 30 only (the allele would add 10 and 20: the reference stays).
    std::string a = gene.substr(0, 150), b = gene.substr(0, 150);
    a[10] = Other(a[10]);
    a[20] = Other(a[20]);
    b[30] = Other(b[30]);
    sam += SamRecord("a", 1, 1, 0, "10=1X9=1X129=", a);
    sam += SamRecord("b", 1, 1, 0, "30=1X119=", b);
    sam += SamRecord("c", 2, 1, 0, "150=", ref.Gene(2, 1).substr(0, 150));
    auto const without = Features(Profile(ref, sam).GetTaxa().at(1));
    for (auto const* name : { "allele_copy_share", "allele_explained_share", "allele_identity_gain" }) EXPECT_EQ(without.at(name), -1.0) << name;
    ref.loader->SetStrainAlleles(table);
    for (size_t threads : { 1, 3 }) {
        SCOPED_TRACE(threads);
        auto const profile = Profile(ref, sam, threads);
        auto const f1 = Features(profile.GetTaxa().at(1));
        auto const f2 = Features(profile.GetTaxa().at(2));
        EXPECT_EQ(f1.at("allele_copy_share"), 1.0);
        EXPECT_NEAR(f1.at("allele_explained_share"), 2.0 / 3, 1e-12);
        EXPECT_NEAR(f1.at("allele_identity_gain"), 2.0 / 300, 1e-12);
        // Taxon 2 has no alleles: nothing explained, as for reads its alleles would not explain (not -1, which would tell
        // a species with other genomes in GTDB from one without, the cluster size); only allele_copy_share says it.
        EXPECT_EQ(f2.at("allele_copy_share"), 0.0);
        EXPECT_EQ(f2.at("allele_explained_share"), 0.0);
        EXPECT_EQ(f2.at("allele_identity_gain"), 0.0);
    }
}

// The polymorphic features. Taxon 1's gene differs from its congener taxon 2's at six ancestry sites (20, 40, ..., 120);
// taxon 1's one allele has the congener's base at 40 (a shared polymorphism) and a third base at 80. Read a has the
// congener's base at 40 (a polymorphic site: a known variant), b at 20 (fixed within the species), c a fourth base at 80
// (polymorphic: no allele's). -1 without the table; a species without alleles (taxon 2) 0, as for the alleles group.
TEST(StrainAlleles, ThePolymorphicFeaturesOfATaxon) {
    std::mt19937 rng(43);
    std::string const gene = test::RandomSequence(600, rng);
    std::string congener = gene;
    for (size_t p : { 20, 40, 60, 80, 100, 120 }) congener[p] = Other(gene[p]);
    auto pick = [](std::string const& taken) {
        for (char const b : std::string("ACGT")) {
            if (taken.find(b) == std::string::npos) return b;
        }
        return 'N';
    };
    char const allele80 = pick({ gene[80], congener[80] });
    char const fourth80 = pick({ gene[80], congener[80], allele80 });
    test::LoadedReference ref({ { 1, { gene } }, { 2, { congener } } }, "polymorphic features");
    ASSERT_EQ(ref.Gene(1, 1), gene);
    species_neighbours::Table neighbours;
    neighbours.Set(1, { { 2, 0.01f } });
    neighbours.Set(2, { { 1, 0.01f } });
    ref.loader->SetSpeciesNeighbours(neighbours);
    auto const table = sa::Table::FromRows({ { 1, 1, { sa::Allele{ 0, 600, { sa::Edit::Substitution(40, sa::BaseCode(congener[40])),
                                                                              sa::Edit::Substitution(80, sa::BaseCode(allele80)) } } } } });
    std::string a = gene.substr(0, 150), b = a, c = a;
    a[40] = congener[40];
    b[20] = congener[20];
    c[80] = fourth80;
    std::string sam = ref.Header();
    sam += SamRecord("a", 1, 1, 0, "40=1X109=", a);
    sam += SamRecord("b", 1, 1, 0, "20=1X129=", b);
    sam += SamRecord("c", 1, 1, 0, "80=1X69=", c);
    sam += SamRecord("d", 2, 1, 0, "150=", congener.substr(0, 150));
    std::vector<char const*> const names{ "polymorphic_known_share", "polymorphic_novel_share", "ancestry_fixed_gain",
                                          "ancestry_fixed_agreement", "allele_sites_per_kb", "ancestry_fixed_share" };
    auto const without = Features(Profile(ref, sam).GetTaxa().at(1));
    for (auto const* name : names) EXPECT_EQ(without.at(name), -1.0) << name;
    // The ancestry sites: 18 visits, 15 with the species' base (a's at 40 and b's at 20 the congener's, c's at 80 neither).
    EXPECT_NEAR(without.at("ancestry_agreement"), 15.0 / 18, 1e-12);
    ref.loader->SetStrainAlleles(table);
    for (size_t threads : { 1, 3 }) {
        SCOPED_TRACE(threads);
        auto const profile = Profile(ref, sam, threads);
        auto const f1 = Features(profile.GetTaxa().at(1));
        // Six polymorphic site visits (40 and 80 by each read): a's at 40 known, c's at 80 novel.
        EXPECT_NEAR(f1.at("polymorphic_known_share"), 1.0 / 6, 1e-12);
        EXPECT_NEAR(f1.at("polymorphic_novel_share"), 1.0 / 6, 1e-12);
        // The fixed sites (20, 60, 100, 120): 11 of 12 the species' base (b's at 20 not), against 15 of 18 at all.
        EXPECT_NEAR(f1.at("ancestry_fixed_gain"), 11.0 / 12 - 15.0 / 18, 1e-12);
        EXPECT_NEAR(f1.at("ancestry_fixed_agreement"), 11.0 / 12, 1e-12);
        EXPECT_NEAR(f1.at("allele_sites_per_kb"), 1000.0 * 6 / 450, 1e-9);
        EXPECT_NEAR(f1.at("ancestry_fixed_share"), 12.0 / 18, 1e-12);
        auto const f2 = Features(profile.GetTaxa().at(2));
        for (auto const* name : names) EXPECT_EQ(f2.at(name), 0.0) << name;
    }
}
