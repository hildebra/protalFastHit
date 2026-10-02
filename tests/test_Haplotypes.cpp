// Unit tests for the phasing of a long-read sample's strains (Profiling/Haplotypes.h): reads clustered by their alleles
// across genes into haplotypes, blocks no read links joined by the strains' shares, errors not taken for strains,
// reads cut between unlikely neighbours, the alleles a read shows (ReadAlleles from VariantHandler::AddAlignment), and
// a strain row called from its own reads (HaplotypeItem).
#include <gtest/gtest.h>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include "Profiling/Haplotypes.h"
#include "Profiling/Strain.h"

using namespace protal;
using haplotypes::ReadRecord;
using haplotypes::Site;

namespace {
    // Sites with two alleles: (gene, position), the reference's base first.
    std::vector<Site> Sites(std::vector<std::pair<uint32_t, uint32_t>> const& at) {
        std::vector<Site> sites;
        for (auto [gene, pos] : at) sites.push_back({ gene, pos, 'A', { 'A', 'G' } });
        return sites;
    }

    // A read's record on `gene` over [start, end) with base `base` at each of `snps` (the reference's A elsewhere).
    ReadRecord Record(uint32_t link, uint32_t gene, uint32_t start, uint32_t end, std::vector<uint32_t> const& snps,
                      uint32_t read_start = 0, char base = 'G') {
        ReadRecord r;
        r.link = link;
        r.gene = gene;
        r.read_start = read_start;
        r.read_end = read_start + (end - start);
        r.alleles.start = start;
        r.alleles.end = end;
        for (auto pos : snps) r.alleles.snps.emplace_back(pos, base);
        return r;
    }

    // `reads` reads of a strain on genes 1 (sites 10, 20, 30) and 2 (sites 5, 15), one record each: the strain with
    // `alt` has the second allele (G) at every site, the other the reference's.
    // These tests phase one or two genes; a run needs three (Settings::min_genes).
    haplotypes::Settings Few() {
        haplotypes::Settings settings;
        settings.min_genes = 1;
        return settings;
    }

    void Strain(std::vector<ReadRecord>& records, uint32_t& link, size_t reads, bool alt, bool genes_linked = true) {
        for (size_t i = 0; i < reads; i++) {
            records.push_back(Record(link, 1, 0, 40, alt ? std::vector<uint32_t>{ 10, 20, 30 } : std::vector<uint32_t>{}, 0));
            if (!genes_linked) link++;
            records.push_back(Record(link, 2, 0, 20, alt ? std::vector<uint32_t>{ 5, 15 } : std::vector<uint32_t>{}, 500));
            link++;
        }
    }
}

TEST(Haplotypes, ReadsAcrossGenesGiveEachStrainItsRow) {
    std::vector<ReadRecord> records;
    uint32_t link = 0;
    Strain(records, link, 14, false);
    Strain(records, link, 6, true);
    auto const sites = Sites({ { 1, 10 }, { 1, 20 }, { 1, 30 }, { 2, 5 }, { 2, 15 } });
    auto const p = haplotypes::Phase(records, sites, 127, {}, Few());
    ASSERT_EQ(p.rows, 2u);
    ASSERT_EQ(p.shares.size(), 2u);
    EXPECT_NEAR(p.shares[0], 0.7, 1e-12);
    EXPECT_NEAR(p.shares[1], 0.3, 1e-12);
    ASSERT_EQ(p.blocks.size(), 1u);  // the reads link both genes
    EXPECT_TRUE(p.blocks[0].phased);
    EXPECT_EQ(p.blocks[0].haplotype_reads, (std::vector<size_t>{ 14, 6 }));
    EXPECT_EQ(p.blocks[0].genes, (std::vector<uint32_t>{ 1, 2 }));
    for (auto const& site : sites) {
        EXPECT_EQ(p.Call(site.gene, site.pos, 0), 'A') << site.gene << ':' << site.pos;
        EXPECT_EQ(p.Call(site.gene, site.pos, 1), 'G') << site.gene << ':' << site.pos;
    }
    EXPECT_EQ(p.Calls(1, 1), 3u);
    EXPECT_EQ(p.Call(1, 11, 1), 0);  // no site
    EXPECT_EQ(p.reads, 20u);
    // Each row's reads: both records of each of its strain's reads.
    ASSERT_EQ(p.row_records.size(), 2u);
    EXPECT_EQ(p.row_records[0].size(), 28u);
    EXPECT_EQ(p.row_records[1].size(), 12u);
    for (auto i : p.row_records[1]) EXPECT_GE(records[i].link, 14u);
}

TEST(Haplotypes, BlocksNoReadLinksJoinByTheStrainsShares) {
    // The genes' records on reads of their own: two blocks, each with the strains at 70:30. The more abundant strain
    // is the first row in both.
    std::vector<ReadRecord> records;
    uint32_t link = 0;
    Strain(records, link, 14, false, false);
    Strain(records, link, 6, true, false);
    auto const sites = Sites({ { 1, 10 }, { 1, 20 }, { 1, 30 }, { 2, 5 }, { 2, 15 } });
    auto p = haplotypes::Phase(records, sites, 127, {}, Few());
    ASSERT_EQ(p.rows, 2u);
    ASSERT_EQ(p.blocks.size(), 2u);
    EXPECT_TRUE(p.blocks[0].phased);
    EXPECT_TRUE(p.blocks[1].phased);
    EXPECT_GT(p.blocks[0].log_odds, 3.0);
    EXPECT_EQ(p.Call(1, 10, 0), 'A');
    EXPECT_EQ(p.Call(2, 5, 0), 'A');
    EXPECT_EQ(p.Call(1, 10, 1), 'G');
    EXPECT_EQ(p.Call(2, 5, 1), 'G');

    // At 50:50 the shares cannot tell which strain of one block is which of the other: the anchor (gene 1's block, of
    // more sites) gives the rows its two strains, gene 2's block keeps the sample's calls.
    records.clear();
    link = 0;
    Strain(records, link, 10, false, false);
    Strain(records, link, 10, true, false);
    p = haplotypes::Phase(records, sites, 127, {}, Few());
    EXPECT_EQ(p.rows, 2u);
    ASSERT_EQ(p.blocks.size(), 2u);
    EXPECT_TRUE(p.blocks[0].phased);
    EXPECT_EQ(p.blocks[0].genes, (std::vector<uint32_t>{ 1 }));
    EXPECT_FALSE(p.blocks[1].phased);
    EXPECT_EQ(p.blocks[1].haplotype_reads, (std::vector<size_t>{ 10, 10 }));
    EXPECT_NEAR(p.blocks[1].log_odds, 0, 1e-9);
    EXPECT_NE(p.Call(1, 10, 0), p.Call(1, 10, 1));
    EXPECT_EQ(p.Call(2, 5, 0), 0);
    EXPECT_EQ(p.Call(2, 5, 1), 0);

    // One block that the reads link is phased within itself, whatever the shares: each row is one strain.
    records.clear();
    link = 0;
    Strain(records, link, 10, false);
    Strain(records, link, 10, true);
    p = haplotypes::Phase(records, sites, 127, {}, Few());
    ASSERT_EQ(p.rows, 2u);
    for (size_t row : { 0, 1 }) {
        char const first = p.Call(1, 10, row);
        ASSERT_TRUE(first == 'A' || first == 'G');
        for (auto const& site : sites) EXPECT_EQ(p.Call(site.gene, site.pos, row), first);
    }
    EXPECT_NE(p.Call(1, 10, 0), p.Call(1, 10, 1));
}

TEST(Haplotypes, ABlockLongerThanItsReadsIsPhasedAlongThem) {
    // 200 sites on one gene, 10 bases apart; reads of two strains at 50:50, each covering 12 consecutive sites at a
    // random place: one block that no read spans. Each row is one strain all along, no halves of two.
    std::mt19937 rng(5);
    std::vector<ReadRecord> records;
    std::vector<std::pair<uint32_t, uint32_t>> at;
    for (uint32_t i = 0; i < 200; i++) at.emplace_back(1, 10 * i);
    for (uint32_t link = 0; link < 400; link++) {
        uint32_t const first = rng() % 189;
        std::vector<uint32_t> snps;
        if (link % 2) {
            for (uint32_t i = first; i < first + 12; i++) snps.push_back(10 * i);
        }
        records.push_back(Record(link, 1, 10 * first, 10 * (first + 12) - 9, snps));
    }
    auto const p = haplotypes::Phase(records, Sites(at), 127, {}, Few());
    ASSERT_EQ(p.blocks.size(), 1u);
    ASSERT_EQ(p.rows, 2u);
    EXPECT_EQ(p.blocks[0].haplotype_reads.size(), 2u);
    for (size_t row : { 0, 1 }) {
        size_t a = 0, g = 0;
        for (auto const& [gene, pos] : at) {
            char const c = p.Call(gene, pos, row);
            a += c == 'A';
            g += c == 'G';
        }
        // A site at an end may lack one strain's reads: the sample's call there.
        EXPECT_TRUE((a >= 190 && g == 0) || (a == 0 && g >= 190)) << "row " << row << ": " << a << " A, " << g << " G";
    }
}

TEST(Haplotypes, TwoStrainsOfTwoGenesAreNoStrainsOfASpecies) {
    // As in ReadsAcrossGenesGiveEachStrainItsRow, but with a run's settings: a second allele on two genes only is a gene
    // from elsewhere, not a strain.
    std::vector<ReadRecord> records;
    uint32_t link = 0;
    Strain(records, link, 14, false);
    Strain(records, link, 6, true);
    auto const p = haplotypes::Phase(records, Sites({ { 1, 10 }, { 1, 20 }, { 1, 30 }, { 2, 5 }, { 2, 15 } }), 127);
    EXPECT_EQ(p.rows, 0u);
    EXPECT_TRUE(p.row_records.empty());
    ASSERT_EQ(p.blocks.size(), 1u);
    EXPECT_FALSE(p.blocks[0].phased);
    EXPECT_EQ(p.blocks[0].haplotype_reads, (std::vector<size_t>{ 14, 6 }));
}

TEST(Haplotypes, ASecondAlleleNoOtherSiteGoesWithIsNoStrain) {
    // One strain; at site 20 four of its 20 reads show G (a systematic error), the other sites agree: no split.
    std::vector<ReadRecord> records;
    for (uint32_t link = 0; link < 20; link++) {
        records.push_back(Record(link, 1, 0, 40, link % 5 == 0 ? std::vector<uint32_t>{ 20 } : std::vector<uint32_t>{}));
    }
    auto const sites = Sites({ { 1, 10 }, { 1, 20 }, { 1, 30 } });
    auto const p = haplotypes::Phase(records, sites, 127, {}, Few());
    EXPECT_EQ(p.rows, 0u);
    ASSERT_EQ(p.blocks.size(), 1u);
    EXPECT_EQ(p.blocks[0].haplotype_reads, (std::vector<size_t>{ 20 }));
}

TEST(Haplotypes, ThreeStrains) {
    // Strains at 50:30:20 on one gene (sites 10, 20, 30, 40): A A A A, G G A A, A A G G.
    std::vector<ReadRecord> records;
    uint32_t link = 0;
    auto add = [&](size_t n, std::vector<uint32_t> const& snps) {
        for (size_t i = 0; i < n; i++) records.push_back(Record(link++, 1, 0, 50, snps));
    };
    add(25, {});
    add(15, { 10, 20 });
    add(10, { 30, 40 });
    auto const sites = Sites({ { 1, 10 }, { 1, 20 }, { 1, 30 }, { 1, 40 } });
    auto const p = haplotypes::Phase(records, sites, 127, {}, Few());
    ASSERT_EQ(p.rows, 3u);
    EXPECT_NEAR(p.shares[0], 0.5, 1e-12);
    EXPECT_NEAR(p.shares[1], 0.3, 1e-12);
    EXPECT_NEAR(p.shares[2], 0.2, 1e-12);
    std::string row0, row1, row2;
    for (uint32_t pos : { 10, 20, 30, 40 }) {
        row0 += p.Call(1, pos, 0);
        row1 += p.Call(1, pos, 1);
        row2 += p.Call(1, pos, 2);
    }
    EXPECT_EQ(row0, "AAAA");
    EXPECT_EQ(row1, "GGAA");
    EXPECT_EQ(row2, "AAGG");
}

TEST(Haplotypes, OwnReadsOnlyAndReadsCutBetweenUnlikelyNeighbours) {
    // A relative's reads (divergence above the MSA's) are left out: here they are all a "second strain".
    std::vector<ReadRecord> records;
    uint32_t link = 0;
    Strain(records, link, 14, false);
    size_t const own = records.size();
    Strain(records, link, 10, true);
    for (size_t i = own; i < records.size(); i++) records[i].divergence = 20;
    auto const sites = Sites({ { 1, 10 }, { 1, 20 }, { 1, 30 }, { 2, 5 }, { 2, 15 } });
    EXPECT_EQ(haplotypes::Phase(records, sites, 10, {}, Few()).rows, 0u);
    EXPECT_EQ(haplotypes::Phase(records, sites, 10, {}, Few()).reads, 14u);
    EXPECT_EQ(haplotypes::Phase(records, sites, 127, {}, Few()).rows, 2u);

    // Every read cut between its two genes: two blocks, as without the links.
    int asked = 0;
    haplotypes::Breaks const always = [&asked](ReadRecord const& a, ReadRecord const& b) {
        asked++;
        return a.gene == 1 && b.gene == 2;
    };
    auto const p = haplotypes::Phase(records, sites, 127, always, Few());
    EXPECT_EQ(p.cut_reads, 24u);
    EXPECT_EQ(p.blocks.size(), 2u);
    EXPECT_EQ(asked, 24);

    // By the gene neighbours: gene 1's 3' end never faces gene 2's 5' end in the taxon's family of 6.
    gene_neighbours::Table table;
    std::istringstream is("# protal gene neighbours: genomes=6 max_gap=3000\n"
                          "clade\tgene\tend\tpartner\tpartner_end\tspecies\tinformative\tgap_median\tgap_min\tgap_max\n"
                          "100\t1\t3\t3\t5\t6\t6\t20\t10\t30\n");
    ASSERT_EQ(table.Read(is), "");
    table.SetLineage(7, { 7, 100 });
    auto const breaks = haplotypes::UnlikelyNeighbours(table, 7);
    ASSERT_TRUE(breaks);
    auto a = Record(0, 1, 0, 40, {}, 0), b = Record(0, 2, 0, 20, {}, 500);
    EXPECT_TRUE(breaks(a, b));    // forward on gene 1, then forward on gene 2: 3' faces 5'
    b.read_start = 5000;          // more than max_gap further on: no neighbours
    EXPECT_FALSE(breaks(a, b));
    b = Record(0, 3, 0, 20, {}, 500);
    EXPECT_FALSE(breaks(a, b));   // gene 3 is the expected neighbour
    EXPECT_FALSE(haplotypes::UnlikelyNeighbours(gene_neighbours::Table{}, 7));
}

TEST(Haplotypes, AReadShowsItsAlleles) {
    // ReadAlleles from AddAlignment: SNPs by their base, N and deleted positions as none, the reference elsewhere.
    std::string const reference = "ACGTTGCAAGGCTTACCGATGACTGAAACCGGTTTACGATCGGTAGCATG";
    VariantHandler handler(reference);
    SamEntry sam;
    sam.m_flag = 0;
    sam.m_pos = 1;
    // 10 matches, a SNP at 10 (T for G), 9 matches, 2 deleted (20, 21), an N at 22 (as M), 27 matches to the end.
    std::string seq = reference.substr(0, 10) + "T" + reference.substr(11, 9) + "N" + reference.substr(23);
    sam.m_seq = seq;
    sam.m_qual = std::string(seq.size(), 'I');
    sam.m_cigar = "10M1X9M2D1M27M";
    ReadAlleles alleles;
    ASSERT_TRUE(handler.AddAlignment(sam, 0, 0, 0, &alleles).has_value());
    EXPECT_EQ(alleles.start, 0u);
    EXPECT_EQ(alleles.end, 50u);
    EXPECT_EQ(alleles.BaseAt(10, 'G'), 'T');
    EXPECT_EQ(alleles.BaseAt(5, reference[5]), reference[5]);
    EXPECT_EQ(alleles.BaseAt(20, reference[20]), 0);
    EXPECT_EQ(alleles.BaseAt(21, reference[21]), 0);
    EXPECT_EQ(alleles.BaseAt(22, reference[22]), 0);
    EXPECT_EQ(alleles.BaseAt(50, 'A'), 0);
}

TEST(Haplotypes, AStrainRowIsCalledFromItsOwnReads) {
    // A sample of two strains of a 10-base gene: strain 1 (6 reads) the reference, strain 2 (4 reads) G at 4 and T at 6;
    // one read of strain 2 has no base at 8 (an N), and none of strain 2 covers 0-1. The sample's own row has R at
    // 4 and K at 6 (both alleles pass); each strain row is its strain's, with a gap where its reads do not reach.
    std::string const reference = "ACGTACGTAC";
    std::vector<ReadRecord> strain1, strain2;
    for (uint32_t i = 0; i < 6; i++) {
        auto r = Record(i, 1, 0, 10, {});
        r.forward = i % 2 == 0;
        strain1.push_back(r);
    }
    for (uint32_t i = 0; i < 4; i++) {
        auto r = Record(10 + i, 1, 2, 10, {});
        r.alleles.snps = { { 4, 'G' }, { 6, 'T' } };
        r.alleles.snp_quals = { 40, 40 };
        r.forward = i % 2 == 0;
        if (i == 0) r.alleles.no_base = { 8 };
        strain2.push_back(r);
    }
    auto pointers = [](std::vector<ReadRecord> const& v) {
        std::vector<ReadRecord const*> out;
        for (auto const& r : v) out.push_back(&r);
        return out;
    };
    auto const one = haplotypes::HaplotypeItem(pointers(strain1), reference, 2, 0.15, 0, 0, false);
    auto const two = haplotypes::HaplotypeItem(pointers(strain2), reference, 2, 0.15, 0, 0, false);
    EXPECT_TRUE(one.first.empty());  // no SNP
    EXPECT_EQ(one.second, CoverageVec(10, 6));
    EXPECT_EQ(two.second, (CoverageVec{ 0, 0, 4, 4, 4, 4, 4, 4, 3, 4 }));
    ASSERT_EQ(two.first.size(), 2u);
    std::vector<ReadRecord> all = strain1;
    all.insert(all.end(), strain2.begin(), strain2.end());
    auto const sample = haplotypes::HaplotypeItem(pointers(all), reference, 2, 0.15, 0, 0, false);
    MSASequenceItems items{ OptionalMSASequenceItem{ sample }, OptionalMSASequenceItem{ one }, OptionalMSASequenceItem{ two } };
    MSAVector msa(3);
    ASSERT_TRUE(MSA(items, reference, msa, 2, 0, std::vector<double>(3, 0.15), false, 0, nullptr, nullptr, 2, 1));
    EXPECT_EQ(std::string(msa[0].begin(), msa[0].end()), "ACGTRCKTAC");
    EXPECT_EQ(std::string(msa[1].begin(), msa[1].end()), "ACGTACGTAC");
    EXPECT_EQ(std::string(msa[2].begin(), msa[2].end()), "--GTGCTTAC");
    // The multi-allelic sites the strains are phased at.
    auto const sites = MultiAllelicSites(sample.first, sample.second, 2, 0, 0.15, false, 0, 4);
    ASSERT_EQ(sites.size(), 2u);
    EXPECT_EQ(sites[0].first, 4u);
    EXPECT_EQ(sites[0].second, (std::vector<char>{ 'A', 'G' }));
}
