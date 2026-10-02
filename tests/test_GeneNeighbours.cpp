// Unit tests for the database's gene neighbours (SequenceUtils/GeneNeighbours.h: gene_neighbours.tsv, which marker
// gene lies next to which in a clade's genomes) and their uses: pairs of mates on two neighbouring genes
// (across_genes::PairAcrossNeighbours, classify::JoinAlignmentPairs), mate guidance past a gene's end
// (GuideMate) and the profiler's adjacent_expected_share and adjacent_unlikely_share.
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
#include "Classify.h"
#include "Profiling/Profiler.h"

using namespace protal;
using gene_neighbours::End;
using gene_neighbours::Verdict;

namespace {
    // Species 1 (family 100, order 200), species 2 (family 101, order 200), species 3 (family 102, order 200),
    // species 5 (family 103, order 200). Family 100: of 6 species, gene 1's 3' end faces gene 2's 5' end in 3 (10-30
    // bases apart, median 20), in the other 3 nothing; gene 2's 5' end, informative in 3 species (fewer than
    // kMinInformative), faces gene 1's 3' end. Family 101: gene 1's 3' end informative in one species, with nothing
    // there. Family 103: gene 1's 3' end faces gene 2's 5' end in 29 of 30, gene 4's in 1. Order 200 (these families
    // and others): gene 1's 3' end faces gene 2's 5' end in 36 of 45, gene 3's in 4, gene 4's in 1, nothing in 4.
    std::string const kTable =
            "# protal gene neighbours: genomes=60 species=45 max_gap=3000 ranks=family,order\n"
            "clade\tgene\tend\tpartner\tpartner_end\tspecies\tinformative\tgap_median\tgap_min\tgap_max\n"
            "100\t1\t3\t2\t5\t3\t6\t20\t10\t30\n"
            "100\t1\t3\t0\t0\t3\t6\t0\t0\t0\n"
            "100\t2\t5\t1\t3\t3\t3\t20\t10\t30\n"
            "101\t1\t3\t0\t0\t1\t1\t0\t0\t0\n"
            "103\t1\t3\t2\t5\t29\t30\t20\t10\t30\n"
            "103\t1\t3\t4\t5\t1\t30\t40\t40\t40\n"
            "200\t1\t3\t2\t5\t36\t45\t20\t10\t30\n"
            "200\t1\t3\t3\t5\t4\t45\t50\t50\t50\n"
            "200\t1\t3\t4\t5\t1\t45\t40\t40\t40\n"
            "200\t1\t3\t0\t0\t4\t45\t0\t0\t0\n";

    gene_neighbours::Table Table() {
        gene_neighbours::Table table;
        std::istringstream is(kTable);
        EXPECT_EQ(table.Read(is), "");
        table.SetLineage(1, { 1, 100, 200, 300 });
        table.SetLineage(2, { 2, 101, 200, 300 });
        table.SetLineage(3, { 3, 102, 200, 300 });
        table.SetLineage(5, { 5, 103, 200, 300 });
        return table;
    }

    std::string Problem(std::string const& line) {
        gene_neighbours::Table table;
        std::istringstream is("clade\tgene\tend\tpartner\tpartner_end\tspecies\tinformative\tgap_median\tgap_min\tgap_max\n" + line);
        return table.Read(is);
    }

    // Taxon 1 with genes 1-3 of 300 random bases each, loaded, with the neighbours of Table().
    struct ThreeGenes {
        std::vector<std::string> genes;
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;

        explicit ThreeGenes(bool with_neighbours = true) {
            static int instances = 0;
            dir = std::filesystem::temp_directory_path() /
                  ("protal_neighbours_" + std::to_string(::getpid()) + "_" + std::to_string(instances++));
            std::mt19937 rng(11);
            std::filesystem::create_directories(dir);
            std::ofstream fna(dir / "reference.fna"), map(dir / "reference.map");
            size_t offset = 0;
            for (int id = 1; id <= 3; id++) {
                std::string seq;
                for (int i = 0; i < 300; i++) seq += "ACGT"[rng() % 4];
                genes.push_back(seq);
                std::string const header = ">1_" + std::to_string(id) + "\n";
                fna << header << seq << '\n';
                map << "1\t" << id << '\t' << offset + header.size() << '\t' << offset + header.size() + seq.size() << '\n';
                offset += header.size() + seq.size() + 1;
            }
            fna.close();
            map.close();
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
            loader->LoadAllGenomes();
            if (with_neighbours) loader->SetGeneNeighbours(Table());
        }
        ~ThreeGenes() {
            loader.reset();
            std::filesystem::remove_all(dir);
        }
    };

    // An ungapped alignment of `length` bases at position `start` of gene `geneid` of taxon 1.
    AlignmentResult Aligned(uint32_t geneid, size_t start, size_t length, bool forward) {
        AlignmentResult ar(0, 1, geneid, static_cast<int32_t>(start), forward);
        auto& info = ar.GetAlignmentInfo();
        info.cigar = std::string(length, 'M');
        info.compressed_cigar = std::to_string(length) + "M";
        info.gene_alignment_start = static_cast<int>(start);
        info.alignment_length = static_cast<uint32_t>(length);
        info.matches = static_cast<uint16_t>(length);
        return ar;
    }

    FastxRecord Record(std::string id, std::string seq) {
        FastxRecord r;
        r.id = std::move(id);
        r.header = "@" + r.id;
        r.sequence = std::move(seq);
        r.quality = std::string(r.sequence.size(), 'I');
        return r;
    }

    // Stands in for the alignment handler: aligns each anchor where it says, ungapped.
    struct UngappedHandler {
        std::vector<CAlignmentAnchor> seen;
        void operator()(AlignmentAnchorList& anchors, AlignmentResultList& results, std::string const& sequence,
                        std::string const&, size_t, std::string&) {
            for (auto const& a : anchors) {
                seen.push_back(a);
                auto const start = static_cast<int32_t>(a.chain.front().genepos) - static_cast<int32_t>(a.chain.front().readpos);
                results.push_back(Aligned(a.geneid, static_cast<size_t>(start), sequence.size(), a.forward));
            }
        }
    };
}

TEST(GeneNeighbours, ReadsTheTableAndItsComment) {
    auto const table = Table();
    EXPECT_EQ(table.Rules(), 10u);
    EXPECT_EQ(table.Clades(), 4u);
    EXPECT_EQ(table.Genomes(), 60u);
    EXPECT_EQ(table.MaxGap(), 3000u);
    EXPECT_EQ(table.BoundSpecies(), 4u);  // clade 300 has no rules, so each lineage is cut to the clades with some
    EXPECT_TRUE(gene_neighbours::Table().Empty());
}

TEST(GeneNeighbours, StopsAtTheFirstBadLine) {
    EXPECT_EQ(Problem("100\t1\t3\t2\t5\t2\t4\t20\t10\t30\n"), "");
    EXPECT_NE(Problem("100\t1\t4\t2\t5\t2\t4\t20\t10\t30\n").find("line 2: end must be 3 or 5"), std::string::npos);
    EXPECT_NE(Problem("100\t1\t3\t2\t0\t2\t4\t20\t10\t30\n").find("partner_end"), std::string::npos);
    EXPECT_NE(Problem("100\t1\t3\t0\t5\t2\t4\t0\t0\t0\n").find("partner_end"), std::string::npos);
    EXPECT_NE(Problem("100\t1\t3\t2\t5\t5\t4\t20\t10\t30\n").find("at most informative"), std::string::npos);
    EXPECT_NE(Problem("100\t1\t3\t2\t5\t2\t4\t5\t10\t30\n").find("gap_min <= gap_median"), std::string::npos);
    EXPECT_NE(Problem("100\t1\t3\t2\t5\t2\t4\t20\t10\n").find("ten"), std::string::npos);
    EXPECT_NE(Problem("100\t1\t3\t2\t5\t2\t4\t20\t10\t30\t7\n").find("ten"), std::string::npos);
    EXPECT_NE(Problem("100\tx\t3\t2\t5\t2\t4\t20\t10\t30\n").find("ten numbers"), std::string::npos);
    EXPECT_NE(Problem("0\t1\t3\t2\t5\t2\t4\t20\t10\t30\n").find("positive"), std::string::npos);
    EXPECT_NE(Problem("100\t1\t3\t2\t5\t2\t4\t20\t10\t30\n100\t1\t3\t2\t5\t1\t4\t20\t10\t30\n").find("listed twice"), std::string::npos);
    EXPECT_NE(Problem("100\t1\t3\t2\t5\t2\t4\t20\t10\t30\n100\t1\t3\t0\t0\t1\t3\t0\t0\t0\n").find("differ in informative"), std::string::npos);
    // Overlapping genes have a negative gap.
    EXPECT_EQ(Problem("100\t1\t3\t2\t5\t2\t4\t-4\t-8\t-1\n"), "");
}

TEST(GeneNeighbours, AFewSpeciesLeanOnTheCladesAbove) {
    auto const table = Table();
    // A clade's share pulled towards its parent's by kPriorSpecies pseudo-species.
    auto smoothed = [](double species, double informative, double parent) {
        return (species + gene_neighbours::kPriorSpecies * parent) / (informative + gene_neighbours::kPriorSpecies);
    };
    double const order2 = 36.0 / 45, order3 = 4.0 / 45, order4 = 1.0 / 45;
    // Species 1: its family saw gene 1's 3' end face gene 2's 5' end in 3 of 6 species, the order in 36 of 45.
    auto expected = table.Assess(1, 1, End::Three, 2, End::Five);
    EXPECT_EQ(expected.verdict, Verdict::Expected);
    EXPECT_EQ(expected.clade, 100u);
    EXPECT_NEAR(expected.share, smoothed(3, 6, order2), 1e-12);  // 0.6
    ASSERT_NE(expected.rule, nullptr);
    EXPECT_EQ(expected.rule->gap_median, 20);
    // ... never gene 3's, which the order has in 4 of 45: unlikely, by the order's gaps.
    auto const never = table.Assess(1, 1, End::Three, 3, End::Five);
    EXPECT_EQ(never.verdict, Verdict::Unlikely);
    EXPECT_EQ(never.clade, 100u);
    EXPECT_NEAR(never.share, smoothed(0, 6, order3), 1e-12);  // 0.03
    ASSERT_NE(never.rule, nullptr);
    EXPECT_EQ(never.rule->gap_median, 50);
    // The other end of gene 2, or another of gene 2's partners: no clade has data on that end.
    EXPECT_EQ(table.Assess(1, 1, End::Five, 2, End::Five).verdict, Verdict::Unknown);
    EXPECT_EQ(table.Assess(1, 1, End::Five, 2, End::Five).clade, 0u);
    // Gene 2's 5' end: informative in 3 species of the family, and no clade above has data. Too few to call a
    // pairing unlikely: gene 1's 3' end there is expected, another partner unknown.
    auto const sparse = table.Assess(1, 2, End::Five, 1, End::Three);
    EXPECT_EQ(sparse.verdict, Verdict::Expected);
    EXPECT_EQ(sparse.clade, 100u);
    EXPECT_DOUBLE_EQ(sparse.share, 1.0);
    EXPECT_EQ(table.Assess(1, 2, End::Five, 3, End::Three).verdict, Verdict::Unknown);
    EXPECT_EQ(table.Assess(1, 2, End::Five, 3, End::Three).clade, 0u);
    // Species 2: its family has the end informative in one species, which has no gene there; the order says
    // most: gene 3 in 4 of 45 is rare for it, though unlikely for species 1, whose family of 6 never has it.
    auto const order = table.Assess(2, 1, End::Three, 3, End::Five);
    EXPECT_EQ(order.verdict, Verdict::Rare);
    EXPECT_EQ(order.clade, 101u);
    EXPECT_NEAR(order.share, smoothed(0, 1, order3), 1e-12);  // 0.067
    EXPECT_NEAR(table.Assess(2, 1, End::Three, 2, End::Five).share, smoothed(0, 1, order2), 1e-12);
    EXPECT_EQ(table.Assess(2, 1, End::Three, 2, End::Five).verdict, Verdict::Expected);  // 0.6
    EXPECT_EQ(table.Assess(2, 1, End::Three, 4, End::Five).verdict, Verdict::Unlikely);  // 0.017
    // Species 3: its family has no data; the order's share is its own.
    EXPECT_EQ(table.Assess(3, 1, End::Three, 2, End::Five).clade, 200u);
    EXPECT_NEAR(table.Assess(3, 1, End::Three, 2, End::Five).share, order2, 1e-12);
    EXPECT_EQ(table.Assess(3, 1, End::Three, 2, End::Three).verdict, Verdict::Unlikely);
    // Species 5: a family of 30 species, one of which has gene 4 there: its own share counts most; unlikely.
    auto const once = table.Assess(5, 1, End::Three, 4, End::Five);
    EXPECT_EQ(once.verdict, Verdict::Unlikely);
    ASSERT_NE(once.rule, nullptr);
    EXPECT_EQ(once.rule->gap_median, 40);
    EXPECT_NEAR(once.share, smoothed(1, 30, order4), 1e-12);  // 0.032
    EXPECT_EQ(table.Assess(5, 1, End::Three, 2, End::Five).verdict, Verdict::Expected);
    // A species of no clade with rules.
    EXPECT_EQ(table.Assess(4, 1, End::Three, 2, End::Five).verdict, Verdict::Unknown);
    // Partners: those of all the species' clades, its family's first, each with the nearest clade's rule and its
    // smoothed share.
    std::vector<gene_neighbours::Partner> partners;
    table.Partners(1, 1, End::Three, partners);
    ASSERT_EQ(partners.size(), 4u);
    std::vector<uint32_t> ids;
    for (auto const& p : partners) ids.push_back(p.rule->partner);
    EXPECT_EQ(ids, (std::vector<uint32_t>{ 0, 2, 3, 4 }));
    EXPECT_EQ(partners[1].rule->informative, 6u);   // the family's rule
    EXPECT_NEAR(partners[1].share, smoothed(3, 6, order2), 1e-12);
    EXPECT_EQ(partners[1].verdict, Verdict::Expected);
    EXPECT_EQ(partners[2].rule->informative, 45u);  // seen by the order only
    EXPECT_EQ(partners[2].verdict, Verdict::Unlikely);
    table.Partners(5, 1, End::Three, partners);
    ASSERT_EQ(partners.size(), 4u);
    EXPECT_NEAR(partners[1].share, smoothed(1, 30, order4), 1e-12);  // gene 4, by the family's rule
    table.Partners(1, 2, End::Five, partners);
    ASSERT_EQ(partners.size(), 1u);
    EXPECT_EQ(partners[0].verdict, Verdict::Expected);
    table.Partners(1, 3, End::Three, partners);
    EXPECT_TRUE(partners.empty());
    // Two genes of a read whose ends may face either way.
    EXPECT_EQ(table.AssessGenes(1, 1, 2), Verdict::Expected);
    EXPECT_EQ(table.AssessGenes(1, 1, 3), Verdict::Unlikely);
    EXPECT_EQ(table.AssessGenes(4, 1, 2), Verdict::Unknown);
}

TEST(GeneNeighbours, WhereAReadLiesOnTheNeighbour) {
    using namespace gene_neighbours;
    EXPECT_EQ(EndAhead(true), End::Three);
    EXPECT_EQ(EndAhead(false), End::Five);
    EXPECT_TRUE(SameStrand(End::Three, End::Five));
    EXPECT_FALSE(SameStrand(End::Three, End::Three));
    EXPECT_TRUE(OrientationOnPartner(true, End::Three, End::Five));    // tandem: the same orientation
    EXPECT_FALSE(OrientationOnPartner(true, End::Three, End::Three));  // convergent: the other
    EXPECT_EQ(ReachToEnd(200, 50, 300, End::Three), 100);  // from the start of a forward mate to the 3' end
    EXPECT_EQ(ReachToEnd(200, 50, 300, End::Five), 250);   // from the end of a reverse mate to the 5' end
    // 0-100 bases past the end, 20 bases apart: the first 80 bases of a gene facing with its 5' end ...
    EXPECT_EQ(PartnerStretch(0, 100, 20, 20, End::Five, 300), (std::pair<int64_t, int64_t>{ 0, 80 }));
    // ... the last 80 of one facing with its 3' end; with a gap range, the widest stretch.
    EXPECT_EQ(PartnerStretch(0, 100, 20, 20, End::Three, 300), (std::pair<int64_t, int64_t>{ 220, 300 }));
    EXPECT_EQ(PartnerStretch(50, 100, 10, 30, End::Five, 300), (std::pair<int64_t, int64_t>{ 20, 90 }));
    // Overlapping genes: the gene begins before the end of the other.
    EXPECT_EQ(PartnerStretch(0, 10, -4, -4, End::Five, 300), (std::pair<int64_t, int64_t>{ 4, 14 }));
    // Nothing reaches the gene, or it ends before.
    auto const none = PartnerStretch(0, 10, 20, 20, End::Five, 300);
    EXPECT_GE(none.first, none.second);
}

TEST(AcrossGenes, MatesOnNeighbouringGenesAreOneFragment) {
    ThreeGenes ref;
    auto& loader = *ref.loader;
    // Mate 1 forward on gene 1 runs to its 3' end, mate 2 reverse on gene 2 to its 5' end: 100 + 10-30 + 150 bases.
    auto const m1 = Aligned(1, 200, 100, true), m2 = Aligned(2, 50, 100, false);
    EXPECT_TRUE(across_genes::PairAcrossNeighbours(m1, m2, loader, 1000));
    EXPECT_FALSE(across_genes::PairAcrossNeighbours(m1, m2, loader, 150));   // longer than the longest fragment
    EXPECT_TRUE(across_genes::PairAcrossNeighbours(m2, m1, loader, 1000));   // mate 2 on gene 1 and mate 1 on gene 2
    // Mate 2 forward on gene 2 runs to its 3' end, which does not face gene 1.
    EXPECT_FALSE(across_genes::PairAcrossNeighbours(m1, Aligned(2, 50, 100, true), loader, 1000));
    // Gene 3 never faces gene 1 in species 1's family.
    EXPECT_FALSE(across_genes::PairAcrossNeighbours(m1, Aligned(3, 50, 100, false), loader, 1000));
    // One gene, or another taxon: not this kind of pair.
    EXPECT_FALSE(across_genes::PairAcrossNeighbours(m1, Aligned(1, 50, 100, false), loader, 1000));
    ThreeGenes plain(false);
    EXPECT_FALSE(across_genes::PairAcrossNeighbours(m1, m2, *plain.loader, 1000));

    // The join pairs them only when asked to.
    for (bool across : { false, true }) {
        AlignmentResultList r1{ m1 }, r2{ m2 };
        PairedAlignmentResultList pairs;
        classify::JoinAlignmentPairs(pairs, r1, r2, loader, false, across);
        ASSERT_EQ(pairs.size(), across ? 1u : 2u);
        EXPECT_TRUE(pairs[0].first.IsSet());
        EXPECT_EQ(pairs[0].second.IsSet(), across);
    }
}

TEST(AcrossGenes, MateGuidanceLooksPastTheGenesEnd) {
    ThreeGenes ref;
    // Mate 1 on the last 50 bases of gene 1, forward; mate 2 on bases 30-80 of gene 2 (reverse, as a mate is): the
    // fragment runs over gene 1's 3' end and the gap into gene 2. Mate 2 has no anchor.
    auto r1 = Record("frag/1", ref.genes[0].substr(250, 50));
    auto r2 = Record("frag/2", KmerUtils::ReverseComplement(ref.genes[1].substr(30, 50)));
    std::string const rev1 = KmerUtils::ReverseComplement(r1.sequence), rev2 = KmerUtils::ReverseComplement(r2.sequence);
    {
        PairedAlignmentResultList pairs{ { Aligned(1, 250, 50, true), AlignmentResult() } };
        UngappedHandler handler;
        MateGuidanceCounts counts;
        EXPECT_TRUE(GuideMate(pairs, {}, {}, r1, r2, rev1, rev2, handler, *ref.loader, counts));
        ASSERT_EQ(handler.seen.size(), 1u);
        EXPECT_EQ(handler.seen[0].geneid, 2u);
        EXPECT_FALSE(handler.seen[0].forward);
        EXPECT_EQ(counts.rescued, 1u);
        EXPECT_EQ(counts.rescued_on_neighbour, 1u);
        // Found across the facing ends: it completes the best candidate.
        ASSERT_EQ(pairs.size(), 1u);
        ASSERT_TRUE(pairs[0].second.IsSet());
        EXPECT_EQ(pairs[0].second.GeneId(), 2u);
        EXPECT_EQ(pairs[0].second.GetAlignmentInfo().gene_alignment_start, 30);
    }
    // Without gene neighbours it is looked for on gene 1 alone, and not found.
    {
        ThreeGenes plain(false);
        auto p1 = Record("frag/1", plain.genes[0].substr(250, 50));
        auto p2 = Record("frag/2", KmerUtils::ReverseComplement(plain.genes[1].substr(30, 50)));
        std::string const prev1 = KmerUtils::ReverseComplement(p1.sequence), prev2 = KmerUtils::ReverseComplement(p2.sequence);
        PairedAlignmentResultList pairs{ { Aligned(1, 250, 50, true), AlignmentResult() } };
        UngappedHandler handler;
        MateGuidanceCounts counts;
        EXPECT_FALSE(GuideMate(pairs, {}, {}, p1, p2, prev1, prev2, handler, *plain.loader, counts));
        EXPECT_EQ(counts.looked_for, 1u);
        EXPECT_EQ(counts.rescued, 0u);
    }
}

TEST(SamReader, KeepsTheHardClipsItDropsFromTheCigar) {
    // Where a long read's record lies on the read: the reader drops the hard clips from the CIGAR and keeps them apart.
    std::string cigar = "120H5S45M30H";
    uint32_t start = 7, end = 7;
    ASSERT_TRUE(NormalizeCigar(cigar, 50, &start, &end));
    EXPECT_EQ(cigar, "5S45M");
    EXPECT_EQ(start, 120u);
    EXPECT_EQ(end, 30u);
    std::string plain = "50M";
    ASSERT_TRUE(NormalizeCigar(plain, 50, &start, &end));
    EXPECT_EQ(start, 0u);
    EXPECT_EQ(end, 0u);
    EXPECT_EQ(profiler::Clip("5S45M", false), 5u);
    EXPECT_EQ(profiler::Clip("120H45M5S30H", true), 35u);
    EXPECT_EQ(profiler::QueryBases("120H5S40M2I3M1D30H"), 45u);
}

TEST(MicrobialProfile, AdjacentGenesOfLinkedReads) {
    ThreeGenes ref;
    auto sam_on = [&](uint32_t gene, int pos, std::string cigar, size_t bases, int flag) {
        SamEntry sam;
        sam.m_qname = "r";
        sam.m_flag = static_cast<FLAG_t>(flag);
        sam.m_rname = "1_" + std::to_string(gene);
        sam.m_pos = pos;
        sam.m_mapq = 60;
        sam.m_cigar = std::move(cigar);
        sam.m_seq = ref.genes[gene - 1].substr(static_cast<size_t>(pos - 1), bases);
        sam.m_qual = std::string(bases, 'I');
        sam.m_alternatives = "*";
        return sam;
    };
    profiler::MicrobialProfile profile(*ref.loader);
    int read = 0;
    // As ProcessMAPQ: every best record is noted for the adjacency, before the filters; the taxon gets its records.
    auto add = [&](SamEntry const& sam, size_t link) {
        int const gene = sam.m_rname.back() - '0';
        profile.NoteLinkedRecord(1, static_cast<uint32_t>(gene), sam, link);
        EXPECT_TRUE(profile.AddSam(1, gene, sam, 1.0, true, read++, true, link));
    };
    // Pair 0: mate 1 forward on gene 1, mate 2 reverse on gene 2: gene 1's 3' end faces gene 2's 5' end.
    add(sam_on(1, 201, "50M", 50, 0x1 | 0x40), 0);
    add(sam_on(2, 31, "50M", 50, 0x1 | 0x80 | 0x10), 0);
    // Pair 1: mate 2 reverse on gene 3, whose 5' end never faces gene 1's 3' end in the family.
    add(sam_on(1, 201, "50M", 50, 0x1 | 0x40), 1);
    add(sam_on(3, 31, "50M", 50, 0x1 | 0x80 | 0x10), 1);
    // Pair 2: both mates on gene 1: not genes next to each other.
    add(sam_on(1, 11, "50M", 50, 0x1 | 0x40), 2);
    add(sam_on(1, 101, "50M", 50, 0x1 | 0x80 | 0x10), 2);
    // A long read: gene 2 (written first, as the best), then gene 1 before it on the read; both forward.
    add(sam_on(2, 1, "120H50M30H", 50, 0x800), 3);
    add(sam_on(1, 251, "50M150H", 50, 0), 3);
    // A long read with gene 1, then gene 2 6 kb further on: more than the table's max_gap apart, so no neighbours.
    add(sam_on(1, 251, "50M6050H", 50, 0), 4);
    add(sam_on(2, 1, "6050H50M", 50, 0x800), 4);
    // A long read whose record of gene 2 the profiler's filters leave out (MAPQ 0): its genes are still neighbours.
    add(sam_on(1, 251, "50M300H", 50, 0), 5);
    auto filtered = sam_on(2, 1, "120H50M180H", 50, 0x800);
    filtered.m_mapq = 0;
    profile.NoteLinkedRecord(1, 2, filtered, 5);
    profile.ApplyRecordEvidence();
    auto const& taxon = profile.GetTaxa().at(1);
    EXPECT_DOUBLE_EQ(taxon.AdjacentExpectedShare(), 3.0 / 4);
    EXPECT_DOUBLE_EQ(taxon.AdjacentUnlikelyShare(), 1.0 / 4);
    std::map<std::string, double> features;
    for (auto const& [name, value] : profiler::TaxonFeatures(taxon)) features[name] = value;
    EXPECT_DOUBLE_EQ(features.at("adjacent_expected_share"), 3.0 / 4);
    EXPECT_DOUBLE_EQ(features.at("adjacent_unlikely_share"), 1.0 / 4);
    // The mean smoothed share of the four links' pairings (gene 2's 0.6 three times, gene 3's 0.03), with one link
    // of 0.5 more.
    double const to2 = (3 + gene_neighbours::kPriorSpecies * 36.0 / 45) / (6 + gene_neighbours::kPriorSpecies);
    double const to3 = (0 + gene_neighbours::kPriorSpecies * 4.0 / 45) / (6 + gene_neighbours::kPriorSpecies);
    EXPECT_NEAR(taxon.AdjacentSupport(), (to2 + to3 + to2 + to2 + 0.5) / 5, 1e-12);
    EXPECT_DOUBLE_EQ(features.at("adjacent_support"), taxon.AdjacentSupport());

    // Without gene neighbours both are 0.
    ThreeGenes plain(false);
    profiler::MicrobialProfile without(*plain.loader);
    SamEntry m1 = sam_on(1, 201, "50M", 50, 0x1 | 0x40), m2 = sam_on(2, 31, "50M", 50, 0x1 | 0x80 | 0x10);
    m1.m_seq = plain.genes[0].substr(200, 50);
    m2.m_seq = plain.genes[1].substr(30, 50);
    without.NoteLinkedRecord(1, 1, m1, 0);
    without.NoteLinkedRecord(1, 2, m2, 0);
    EXPECT_TRUE(without.AddSam(1, 1, m1, 1.0, true, 0, true, 0));
    EXPECT_TRUE(without.AddSam(1, 2, m2, 1.0, true, 1, true, 0));
    without.ApplyRecordEvidence();
    EXPECT_DOUBLE_EQ(without.GetTaxa().at(1).AdjacentExpectedShare(), 0);
    EXPECT_DOUBLE_EQ(without.GetTaxa().at(1).AdjacentUnlikelyShare(), 0);
    EXPECT_DOUBLE_EQ(without.GetTaxa().at(1).AdjacentSupport(), 0.5);
}
