// Unit tests of ColumnWeights.h and ColumnWeightsBuild.h: the weight scale, the alignment maps and runs, a hand-made
// family's columns, the table's round trip and expansion, what a record's mismatches weigh (with the codon
// classification), the "weights" features through the profiler, and the ancestry sites polarised by the family's
// consensus.
#include <gtest/gtest.h>
#include <cmath>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "Core/AlignmentStrategy.h"
#include "Profiling/Profiler.h"
#include "SequenceUtils/AncestrySites.h"
#include "SequenceUtils/ColumnWeights.h"
#include "SequenceUtils/ColumnWeightsBuild.h"
#include "TestReference.h"
#include "TestUtil.h"

using namespace protal;
namespace cw = protal::column_weights;

namespace {
    // A coding sequence of `codons` codons without stops (the frame matters for the amino-acid votes).
    std::string CodingSequence(size_t codons, std::mt19937& rng) {
        static constexpr char kBases[] = "ACGT";
        std::string seq;
        while (seq.size() < 3 * codons) {
            std::string codon{ kBases[rng() % 4], kBases[rng() % 4], kBases[rng() % 4] };
            if (cw::AminoAcid(cw::BaseCode(codon[0]), cw::BaseCode(codon[1]), cw::BaseCode(codon[2])) == '*') continue;
            seq += codon;
        }
        return seq;
    }

    char Other(char c, size_t k = 0) {
        std::string others;
        for (char const b : std::string("ACGT")) {
            if (b != c) others += b;
        }
        return others[k % 3];
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
        auto const path = ref.Write("weights" + std::to_string(threads) + ".sam", sam);
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
}

// The scale: -log(1 - p) in half nats with the pseudocounts (k + 1)/(n + 2), 0 without an estimate, capped at 15.
TEST(ColumnWeights, TheWeightIsTheLogOddsOfTheConservation) {
    EXPECT_EQ(cw::CodeOf(0, 0), cw::kNoCode);
    EXPECT_EQ(cw::CodeOf(1, 1), 2u);     // 2/3: 1.10 nats
    EXPECT_EQ(cw::CodeOf(9, 9), 5u);     // 10/11: 2.40 nats
    EXPECT_EQ(cw::CodeOf(99, 99), 9u);   // 100/101: 4.62 nats: a conserved column
    EXPECT_EQ(cw::CodeOf(999, 999), 14u);
    EXPECT_EQ(cw::CodeOf(100000, 100000), cw::kMaxCode);
    EXPECT_EQ(cw::CodeOf(1, 9), 1u);     // 2/11: 0.2 nats, rounded to the least code
    EXPECT_EQ(cw::CodeOfShare(0.99, 1), 9u);
    EXPECT_DOUBLE_EQ(cw::Nats(9), 4.5);
    EXPECT_GE(cw::kConservedCode, cw::CodeOf(90, 91));  // 90 of 91 agreeing (98.9%) is just below
    EXPECT_EQ(cw::AminoAcid(0, 3, 2), 'M');             // ATG
    EXPECT_EQ(cw::AminoAcid(3, 0, 0), '*');             // TAA
    EXPECT_EQ(cw::AminoAcid(4, 0, 0), 'X');
}

// MapTo: per reference position the copy's position, through substitutions and an indel; the runs of such a map.
TEST(ColumnWeights, TheMapFollowsTheAlignment) {
    std::mt19937 rng(21);
    auto const reference = CodingSequence(200, rng);
    std::string copy = reference;
    copy[100] = Other(copy[100]);
    copy.erase(300, 3);       // the copy lacks the reference's 300-302
    copy.insert(450, "GGG");  // and has three bases the reference lacks before 450
    WFA2Wrapper2 aligner(4, 6, 2, 0);
    std::string cigar;
    auto const map = cw::MapTo(reference, copy, 0.2, aligner, cigar);
    ASSERT_EQ(map.size(), reference.size());
    EXPECT_EQ(map[50], 50);
    EXPECT_EQ(map[100], 100);
    EXPECT_EQ(map[299], 299);
    EXPECT_EQ(map[303], 300);   // after the deletion the copy runs three behind
    EXPECT_EQ(map[500], 500);   // and the insertion brings it back
    int gaps = 0;
    for (size_t i = 300; i < 303; i++) gaps += map[i] < 0;
    EXPECT_EQ(gaps, 3);
    // The copy's columns, as the table stores them: the copy positions in runs.
    std::vector<int32_t> copy_to_column(copy.size(), -1);
    for (size_t r = 0; r < map.size(); r++) {
        if (map[r] >= 0) copy_to_column[static_cast<size_t>(map[r])] = static_cast<int32_t>(r);
    }
    auto const runs = cw::RunsOf(copy_to_column);
    ASSERT_GE(runs.size(), 3u);
    EXPECT_EQ(runs[0], (cw::Run{ 0, 0, 300 }));
    EXPECT_EQ(runs[1].begin, 300u);
    EXPECT_EQ(runs[1].column, 303u);
    // Identical copies map one to one without an alignment; a copy too far away maps nowhere.
    EXPECT_EQ(cw::MapTo(reference, reference, 0.2, aligner, cigar)[123], 123);
    std::string far = reference;
    for (size_t i = 0; i < far.size(); i += 3) far[i] = Other(far[i]);
    EXPECT_TRUE(cw::MapTo(reference, far, 0.2, aligner, cigar).empty());
}

// A family of two genera: the columns where the species of a genus vary get a low within weight, those where the two
// genus references differ a low among weight, the rest high; the amino-acid weight follows the codons; every copy maps.
TEST(ColumnWeights, AFamilysColumnsAreCountedFromItsCopies) {
    std::mt19937 rng(7);
    auto const base = CodingSequence(200, rng);  // 600 bases
    // Genus 10: species 1 (its reference by hash, whichever), 2, 3 vary at 100 (every species another base) and at 101
    // (a synonymous change at a third position would still be a base change: the aa weight sees the codon); genus 20:
    // species 4 and 5, which differ from genus 10 at 200 and 201 (a genus-level difference) and at 300 (species 5 only).
    std::vector<uint32_t> taxids{ 1, 2, 3, 4, 5 };
    std::vector<uint32_t> genus{ 10, 10, 10, 20, 20 };
    std::vector<std::string> seqs(5, base);
    seqs[1][100] = Other(base[100], 0);
    seqs[2][100] = Other(base[100], 1);
    seqs[3][200] = Other(base[200], 0);
    seqs[4][200] = Other(base[200], 0);
    seqs[3][201] = Other(base[201], 0);
    seqs[4][201] = Other(base[201], 0);
    seqs[4][300] = Other(base[300], 0);
    cw::Settings settings;
    cw::Stats stats;
    std::vector<cw::Family> families;
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<cw::Run>>> copies;
    WFA2Wrapper2 aligner(4, 6, 2, 0);
    ASSERT_TRUE(cw::BuildFamily(99, 7, taxids, seqs, genus, settings, aligner, stats, families, copies));
    ASSERT_EQ(families.size(), 1u);
    auto const& f = families[0];
    EXPECT_EQ(f.family, 99u);
    EXPECT_EQ(f.gene, 7u);
    EXPECT_EQ(f.genera, 2u);
    EXPECT_EQ(f.Length(), base.size());
    EXPECT_EQ(copies.size(), 5u);
    EXPECT_EQ(stats.alignments, 4u);  // one genus reference and three species against theirs
    EXPECT_EQ(stats.unaligned, 0u);
    // The reference's columns are the family reference's; the family reference is one of the five (a base sequence
    // apart from its own substitutions), so column 50 is conserved everywhere.
    EXPECT_GT(f.within[50], f.within[100]) << "species vary at 100";
    EXPECT_GT(f.among[50], f.among[200]) << "the genera differ at 200";
    EXPECT_LE(f.among[100], f.among[50]) << "the genus 10 reference may be a species that varies at 100";
    // Only species 5 differs at 300: the genus references agree there unless species 5 is genus 20's reference (the
    // member with the least hash, as BuildFamily chooses it).
    bool const five_is_reference = cw::detail::Mix(5 ^ cw::detail::Mix(20)) < cw::detail::Mix(4 ^ cw::detail::Mix(20));
    if (five_is_reference) EXPECT_LT(f.among[300], f.among[50]);
    else EXPECT_EQ(f.among[300], f.among[50]);
    EXPECT_EQ(f.reference == 4 || f.reference == 5, five_is_reference ? f.reference == 5 : f.reference == 4);
    EXPECT_NE(f.within[100], cw::kNoCode);
    // Two genus references are no outgroup: no consensus base (a majority of at least two among three or more).
    EXPECT_EQ(f.consensus[50], cw::kNoBase);
    EXPECT_EQ(f.consensus[200], cw::kNoBase);
    EXPECT_EQ(f.aa.size(), base.size() / 3);
    EXPECT_NE(f.aa[50 / 3], cw::kNoCode);
    // Every copy maps onto all its positions in one run (substitutions only).
    for (auto const& [taxid, gene, row, runs] : copies) {
        EXPECT_EQ(gene, 7u);
        EXPECT_EQ(row, 0u);
        ASSERT_EQ(runs.size(), 1u) << taxid;
        EXPECT_EQ(runs[0], (cw::Run{ 0, 0, static_cast<uint16_t>(base.size()) })) << taxid;
    }
    // The same through ScanGene with a family id per copy, and a copy without a family left out.
    std::vector<cw::Family> families2;
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<cw::Run>>> copies2;
    cw::Stats stats2;
    std::vector<uint32_t> family{ 99, 99, 99, 99, 0 };
    cw::ScanGene(7, taxids, seqs, genus, family, settings, 2, stats2, families2, copies2);
    ASSERT_EQ(families2.size(), 1u);
    EXPECT_EQ(copies2.size(), 4u);
    EXPECT_EQ(stats2.genes, 1u);
    // A third genus (species 6 and 7, differing from the base at 200 like genus 20, and at 400 on their own) gives a
    // consensus: the base's base where two of three genus references agree, none where all three differ.
    std::vector<uint32_t> taxids3 = taxids, genus3 = genus;
    std::vector<std::string> seqs3 = seqs;
    for (uint32_t t : { 6u, 7u }) {
        taxids3.push_back(t);
        genus3.push_back(30);
        seqs3.push_back(base);
        seqs3.back()[200] = Other(base[200], 0);
        seqs3.back()[201] = Other(base[201], 0);
        seqs3.back()[400] = Other(base[400], 1);
    }
    seqs3[0][500] = Other(base[500], 0);  // genus 10 alone differs at 500, whichever its reference
    seqs3[1][500] = Other(base[500], 0);
    seqs3[2][500] = Other(base[500], 0);
    std::vector<cw::Family> families3;
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<cw::Run>>> copies3;
    cw::Stats stats3;
    ASSERT_TRUE(cw::BuildFamily(99, 7, taxids3, seqs3, genus3, settings, aligner, stats3, families3, copies3));
    auto const& g = families3[0];
    EXPECT_EQ(g.genera, 3u);
    EXPECT_EQ(g.consensus[50], cw::BaseCode(base[50]));
    EXPECT_EQ(g.consensus[200], cw::BaseCode(Other(base[200], 0)));  // genera 20 and 30 against 10
    EXPECT_EQ(g.consensus[500], cw::BaseCode(base[500]));             // genera 20 and 30 keep the base's base
    EXPECT_EQ(g.consensus[400], cw::BaseCode(base[400]));             // genus 30 alone differs
}

// The table round-trips through its text, finds copies and expands their columns.
TEST(ColumnWeights, TheTableRoundTripsAndExpands) {
    cw::Family f;
    f.family = 5;
    f.gene = 2;
    f.reference = 11;
    f.genera = 3;
    f.within = { 3, 9, 15, 0, 7, 7 };
    f.among = { 2, 2, 12, 0, 1, 15 };
    f.aa = { 8, 4 };
    f.consensus = { 0, 1, 2, 4, 3, 0 };
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t, std::vector<cw::Run>>> copies{
        { 11, 2, 0, { cw::Run{ 0, 0, 6 } } },
        { 12, 2, 0, { cw::Run{ 0, 1, 2 }, cw::Run{ 3, 3, 3 } } },  // position 2 of copy 12 maps nowhere
        { 13, 2, 0, { cw::Run{ 0, 0, 6 } } },
    };
    auto const table = cw::Table::FromRows({ f }, copies);
    EXPECT_EQ(table.Families(), 1u);
    EXPECT_EQ(table.Copies(), 3u);
    EXPECT_EQ(table.Species(), 3u);
    EXPECT_EQ(table.ColumnCount(), 6u);
    std::ostringstream os;
    table.Write(os);
    cw::Table back;
    std::istringstream is(os.str());
    EXPECT_EQ(back.Read(is), "") << os.str();
    EXPECT_EQ(back.Copies(), 3u);
    auto const view = back.Find(12, 2);
    ASSERT_TRUE(view.has_value());
    EXPECT_EQ(view->family->within, f.within);
    EXPECT_EQ(view->family->consensus, f.consensus);
    EXPECT_EQ(view->n_runs, 2u);
    EXPECT_FALSE(back.Find(14, 2).has_value());
    EXPECT_FALSE(back.Find(12, 4).has_value());
    auto const columns = back.Expand(12, 2, 6);
    ASSERT_TRUE(columns);
    EXPECT_EQ(columns->within, (std::vector<uint8_t>{ 9, 15, 0, 0, 7, 7 }));
    EXPECT_EQ(columns->among, (std::vector<uint8_t>{ 2, 12, 0, 0, 1, 15 }));
    EXPECT_EQ(columns->aa, (std::vector<uint8_t>{ 8, 8, 0, 4, 4, 4 }));
    EXPECT_EQ(columns->consensus, (std::vector<uint8_t>{ 1, 2, 4, 4, 3, 0 }));
    EXPECT_EQ(columns->Known(), 4u);
    EXPECT_FALSE(back.Expand(14, 2, 6));
    // Malformed text is refused with the line: a copy of a gene its family row does not cover, or no F line first.
    cw::Table bad;
    std::istringstream bad_is("C\t1\t2\t5\t0:0:6\n");
    EXPECT_NE(bad.Read(bad_is).find("line 1"), std::string::npos);
    std::istringstream other_gene(os.str() + "C\t12\t3\t5\t0:0:6\n");
    EXPECT_NE(bad.Read(other_gene).find("line 6"), std::string::npos);
}

// What a record's mismatches weigh: the columns with an estimate, the mismatches' weights, the conserved ones, and
// the codon classification against the copy.
TEST(ColumnWeights, ARecordsMismatchesAreWeighedByTheirColumns) {
    std::string const copy = "ATGAAACCCGGGTTTATG";  // M K P G F M
    cw::Columns columns;
    columns.within = { 9, 9, 9, 2, 2, 2, 0, 0, 0, 12, 12, 12, 1, 1, 1, 9, 9, 9 };
    columns.among = { 5, 5, 5, 1, 1, 1, 0, 0, 0, 3, 3, 3, 1, 1, 1, 5, 5, 5 };
    columns.aa = { 9, 9, 9, 2, 2, 2, 0, 0, 0, 12, 12, 12, 1, 1, 1, 9, 9, 9 };
    columns.consensus.assign(18, 0);
    // The read: the first codon's A->C (ATG -> CTG: M -> L, at a conserved codon), the second codon's third base A->G
    // (AAA -> AAG: synonymous), the third codon's C->T at a column without an estimate (not counted), the fourth's
    // G->A (GGG -> AGG: G -> R, at a conserved column 9 of within 12 and aa 12).
    std::string read = copy;
    read[0] = 'C';
    read[5] = 'G';
    read[6] = 'T';
    read[9] = 'A';
    auto const c = cw::Count(columns, "1X4=1X1X2=1X8=", 1, read, copy);
    EXPECT_EQ(c.columns, 18u);
    EXPECT_EQ(c.aligned, 15u);
    EXPECT_EQ(c.within_aligned, 3u * 9 + 3 * 2 + 3 * 12 + 3 * 1 + 3 * 9);
    EXPECT_EQ(c.mismatches, 3u);
    EXPECT_EQ(c.within_mismatches, 9u + 2 + 12);
    EXPECT_EQ(c.among_mismatches, 5u + 1 + 3);
    EXPECT_EQ(c.conserved_mismatches, 2u);
    EXPECT_EQ(c.synonymous, 1u);
    EXPECT_EQ(c.nonsynonymous, 2u);
    EXPECT_EQ(c.nonsynonymous_conserved, 2u);
    // A record without a sequence counts nothing; one past the copy's end stops at it; a soft clip and a deletion walk.
    EXPECT_EQ(cw::Count(columns, "18M", 1, "*", copy).aligned, 0u);
    EXPECT_EQ(cw::Count(columns, "3S15M", 4, "GGG" + copy.substr(3), copy).aligned, 12u);
    auto const d = cw::Count(columns, "3=3D12=", 1, copy.substr(0, 3) + copy.substr(6), copy);
    EXPECT_EQ(d.aligned, 12u);
    EXPECT_EQ(d.mismatches, 0u);
}

// The features through the profiler: -1 without the table; with it the ratio, the rates and the weighted agreement.
TEST(ColumnWeights, TheFeaturesComeFromTheRecords) {
    std::mt19937 rng(9);
    auto const gene = CodingSequence(200, rng);  // 600 bases
    test::LoadedReference ref({ { 1, { gene } }, { 2, { gene } } });  // species 2 a congener with the same copy
    cw::Family f;
    f.family = 5;
    f.gene = 1;
    f.reference = 1;
    f.genera = 2;
    f.within.assign(600, 9);
    f.among.assign(600, 9);
    f.aa.assign(200, 9);
    f.consensus.resize(600);
    for (size_t i = 0; i < 600; i++) f.consensus[i] = cw::BaseCode(gene[i]);
    for (size_t i = 100; i < 110; i++) f.within[i] = 1;  // a variable stretch
    auto const table = cw::Table::FromRows({ f }, { { 1, 1, 0, { cw::Run{ 0, 0, 600 } } } });
    // Reads a and b carry one mismatch each: a's at a variable column (103), b's at a conserved one (250).
    std::string a = gene.substr(0, 150), b = gene.substr(200, 150);
    a[103] = Other(a[103]);
    b[50] = Other(b[50]);
    std::string sam = ref.Header();
    sam += SamRecord("a", 1, 1, 0, "103=1X46=", a);
    sam += SamRecord("b", 1, 1, 200, "50=1X99=", b);
    auto const without = Features(Profile(ref, sam).GetTaxa().at(1));
    for (auto const* name : { "conserved_mismatch_ratio", "conserved_mismatch_rate", "ancestry_agreement_weighted",
                              "nonsynonymous_share", "nonsynonymous_conserved_rate", "column_weight_coverage" }) {
        EXPECT_EQ(without.at(name), -1.0) << name;
    }
    ref.loader->SetColumnWeights(table);
    for (size_t threads : { 1, 3 }) {
        SCOPED_TRACE(threads);
        auto const features = Features(Profile(ref, sam, threads).GetTaxa().at(1));
        // Mean within of the mismatches (1 + 9)/2 = 5 over the mean of the aligned columns ((290 x 9 + 10 x 1)/300).
        double const aligned_mean = (290.0 * 9 + 10.0 * 1) / 300.0;
        EXPECT_NEAR(features.at("conserved_mismatch_ratio"), 5.0 / aligned_mean, 1e-9);
        EXPECT_NEAR(features.at("conserved_mismatch_rate"), 1000.0 * 1 / 300, 1e-9);
        EXPECT_DOUBLE_EQ(features.at("column_weight_coverage"), 1.0);
        EXPECT_GE(features.at("nonsynonymous_share"), 0.0);
        EXPECT_LE(features.at("nonsynonymous_share"), 1.0);
        EXPECT_EQ(features.at("ancestry_agreement_weighted"), -1.0);  // no congener differs: no site
    }
    // A taxon whose copy has no row in the table: the features known (the table exists) but empty.
    std::string const other = ref.Header() + SamRecord("c", 2, 1, 0, "150=", gene.substr(0, 150));
    auto const empty = Features(Profile(ref, other).GetTaxa().at(2));
    EXPECT_EQ(empty.at("conserved_mismatch_ratio"), 0.0);
    EXPECT_EQ(empty.at("column_weight_coverage"), 0.0);
}

// The weighted site shift in the alignment: a read that fits two species alike by the references alone (four
// differences each) is settled for the species whose differences lie at columns the family's genera change freely
// (among code 1: each half a difference off), against the one whose differences lie at conserved columns (code 9).
// Without the table neither candidate shifts; a copy without strain alleles gets the discount alone.
TEST(ColumnWeights, TheShiftDiscountsDifferencesAtVariableColumns) {
    std::mt19937 rng(41);
    std::string const base = CodingSequence(200, rng);  // 600 bases
    auto substituted = [&base](std::vector<size_t> const& positions) {
        std::string s = base;
        for (size_t const p : positions) s[p] = Other(s[p]);
        return s;
    };
    std::string const rep1 = substituted({ 105, 125, 145, 165 }), rep2 = substituted({ 110, 130, 150, 170 });
    test::LoadedReference ref({ { 1, { rep1 } }, { 2, { rep2 } } }, "weighted shift");
    cw::Family f1, f2;
    f1.family = 10;
    f2.family = 20;
    f1.gene = f2.gene = 1;
    f1.reference = 1;
    f2.reference = 2;
    f1.genera = f2.genera = 5;
    for (auto* f : { &f1, &f2 }) {
        f->within.assign(600, 9);
        f->among.assign(600, 9);
        f->aa.assign(200, 9);
        f->consensus.assign(600, cw::kNoBase);
    }
    for (size_t const p : { 105, 125, 145, 165 }) f1.among[p] = 1;  // taxon 1 differs from the read where the family drifts
    auto const table = cw::Table::FromRows({ f1, f2 }, { { 1, 1, 0, { cw::Run{ 0, 0, 600 } } }, { 2, 1, 1, { cw::Run{ 0, 0, 600 } } } });
    std::string const read = base.substr(100, 150);
    std::string const rev = KmerUtils::ReverseComplement(read);
    WFA2Wrapper2 aligner(4, 6, 2, 1000);
    std::string id = "r";
    for (bool weights : { false, true }) {
        SCOPED_TRACE(weights);
        if (weights) ref.loader->SetColumnWeights(table);
        SimpleAlignmentHandler handler(*ref.loader, aligner, 31, 3, 0.9, false);
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
        handler(anchors, results, read, rev, 3, id);
        ASSERT_EQ(results.size(), 2u);
        auto const& on1 = results.front().Taxid() == 1 ? results.front() : results.back();
        auto const& on2 = results.front().Taxid() == 2 ? results.front() : results.back();
        EXPECT_EQ(on1.GetAlignmentInfo().mismatches, 4);
        EXPECT_EQ(on2.GetAlignmentInfo().mismatches, 4);
        EXPECT_EQ(handler.m_unsure_reads, 1u);  // a tie by the references: unsure either way
        if (!weights) {
            EXPECT_EQ(on1.GetAlignmentInfo().site_shift, 0);
            EXPECT_EQ(on2.GetAlignmentInfo().site_shift, 0);
            EXPECT_EQ(handler.m_settled_reads, 0u);
            continue;
        }
        EXPECT_EQ(on1.GetAlignmentInfo().site_shift, -4);  // four half differences off
        EXPECT_EQ(on2.GetAlignmentInfo().site_shift, 0);   // conserved columns: nothing off
        EXPECT_EQ(results.front().Taxid(), 1u);
        EXPECT_GT(on1.AlignmentScore(), on2.AlignmentScore());
        EXPECT_EQ(handler.m_allele_scored, 0u);  // no strain alleles: the weights alone
        EXPECT_EQ(handler.m_allele_shifted, 1u);
    }
}

// The family's consensus polarises the sites of a species with fewer than three congeners: where the congener's base
// is the family's, the species' is derived (a site); where the congener's own base is new, it is not.
TEST(ColumnWeights, TheFamilyConsensusPolarisesASmallGenus) {
    std::mt19937 rng(31);
    auto const own = CodingSequence(300, rng);
    std::string congener = own;
    congener[100] = Other(own[100], 0);  // the family carries the congener's base: the species changed here
    congener[200] = Other(own[200], 0);  // the family carries the species' base: the congener changed here
    congener[300] = Other(own[300], 0);  // the family carries a third base: unknown
    auto comparison = ancestry::CompareCopies(own, congener);
    ASSERT_TRUE(comparison.has_value());
    comparison->sites.congener = 2;
    std::vector<ancestry::Comparison> comparisons{ std::move(*comparison) };
    std::vector<uint8_t> outgroup(own.size());
    for (size_t i = 0; i < own.size(); i++) outgroup[i] = cw::BaseCode(own[i]);
    outgroup[100] = cw::BaseCode(congener[100]);
    outgroup[300] = cw::BaseCode(Other(own[300], 1));
    outgroup[400] = cw::kNoBase;
    auto const plain = ancestry::Consensus(own.size(), comparisons);
    EXPECT_EQ(plain.positions, (std::vector<uint16_t>{ 100, 200, 300 }));
    auto const polarised = ancestry::Consensus(own.size(), comparisons, &outgroup);
    EXPECT_EQ(polarised.positions, (std::vector<uint16_t>{ 100 }));
    EXPECT_EQ(polarised.bases, (std::vector<uint8_t>{ cw::BaseCode(congener[100]) }));
    // Without a consensus base at a position (kNoBase) the fallback stands.
    std::string other = own;
    other[400] = Other(own[400], 0);
    auto comparison2 = ancestry::CompareCopies(own, other);
    ASSERT_TRUE(comparison2.has_value());
    comparison2->sites.congener = 3;
    std::vector<ancestry::Comparison> comparisons2{ std::move(*comparison2) };
    EXPECT_EQ(ancestry::Consensus(own.size(), comparisons2, &outgroup).positions, (std::vector<uint16_t>{ 400 }));
}
