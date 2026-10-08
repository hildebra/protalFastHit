// Unit tests for the strain alleles (2026-10-08, docs/claude/2026-10-08-strain-alleles): the genomes that may give them
// (the accession hash, the same in Python), the table and its file, an allele's edits from its alignment against the
// representative, the build's sample and choice, what an allele explains of a read, the alignment score's shift, the
// handler's scores, and the profiler's alleles features.
#include <gtest/gtest.h>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include "Core/AlignmentStrategy.h"
#include "Profiling/Profiler.h"
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

// The build's sample of a copy's alleles (the least hashes) and its choice (farthest first) do not depend on the order
// the alleles come in; an allele as far as the nearest congener's copy, of a genome outside the share, identical to the
// representative or seen before is not offered.
TEST(StrainAlleles, TheBuildsSampleAndChoiceDoNotDependOnTheOrder) {
    using E = sa::Edit;
    sa::Allele const a{ 0, 100, { E::Substitution(10, 0), E::Substitution(20, 1), E::Substitution(30, 2) } };  // 3 from the rep
    sa::Allele const b{ 0, 100, { E::Substitution(10, 0) } };                                                   // 1, near a
    sa::Allele const c{ 0, 100, { E::Substitution(40, 3), E::Substitution(50, 3) } };                           // 2, far from a
    EXPECT_EQ(sa::Select({ { 1, a }, { 2, b }, { 3, c } }, 2), (std::vector<sa::Allele>{ a, c }));
    EXPECT_EQ(sa::Select({ { 1, a }, { 2, b }, { 3, c } }, 5), (std::vector<sa::Allele>{ a, c, b }));
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

TEST(StrainAlleles, TheScoreCountsTheShiftAsMatches) {
    AlignmentInfo info;
    info.matches = 97;
    info.mismatches = 3;
    EXPECT_EQ(info.Score(), -12);
    info.allele_shift = -2;
    EXPECT_EQ(info.Score(), -4);
    EXPECT_EQ(info.Score(2, 3, 1, 2), 97 * 2 - 3 * 3 + 2 * 5);
    info.Reset();
    EXPECT_EQ(info.allele_shift, 0);
}

// A read of a known strain of taxon 1 (its allele: the gene before taxon 1's own six differences) fits taxon 2 better by
// the references alone (four differences against six); with the alleles it scores on taxon 1 as on the strain's gene.
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
        EXPECT_EQ(on1.GetAlignmentInfo().allele_shift, alleles ? -6 : 0);
        EXPECT_EQ(handler.m_allele_scored, alleles ? 1u : 0u);
        EXPECT_EQ(handler.m_allele_shifted, alleles ? 1u : 0u);
    }
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
