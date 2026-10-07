// Unit tests for the gene windows of the alignment stage: the anchor finder decodes only the stretch of a
// gene along an anchor's links (ChainAnchorFinder::GeneAround) and extends the anchor there, and must
// get what it gets from the whole gene, including for reads that run past either end of the gene, chains
// with an indel between their links, genes shorter than the read and genes on the heap (> 4096 bases).
// In an AddressSanitizer build a read outside the window is an error (GeneSequence poisons it).
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <memory>
#include <random>
#include <string>
#include <vector>
#include <unistd.h>
#include "SequenceUtils/GenomeLoader.h"
#include "Hash/KmerLookup.h"
#include "Core/AlignmentStrategy.h"
#include "Core/ChainAnchorFinder.h"
#include "TestReference.h"

using namespace protal;

namespace {
    // Random genes of taxid 1 of the given lengths (gene ids 1..), loaded.
    struct Reference : test::LoadedReference {
        std::vector<std::string> const& genes = LoadedReference::genes.at(1);

        explicit Reference(std::vector<size_t> const& lengths) : LoadedReference({ { 1, Genes(lengths) } }, "windows") {}

        static std::vector<std::string> Genes(std::vector<size_t> const& lengths) {
            std::mt19937 rng(31);
            return test::RandomGenes(lengths, rng);
        }
    };

    using Finder = ChainAnchorFinder<KmerLookupSM>;

    // An anchor of gene `gene` with these links, extended by the finder's own window and, for the
    // reference, by the whole gene decoded: the two must be the same.
    void ExpectSameExtension(Finder& finder, Reference& ref, uint32_t gene, std::string const& read, ChainList chain, char const* what) {
        ChainAlignmentAnchor windowed(1, gene, true), whole(1, gene, true);
        windowed.chain = chain;
        whole.chain = chain;
        finder.ExtendAnchor(windowed, read);
        auto const sequence = ref.loader->GetGenome(1).GetGeneOMP(gene).Sequence();
        finder.ExtendAnchor(whole, read, sequence);
        ASSERT_EQ(windowed.chain.size(), whole.chain.size()) << what;
        for (size_t i = 0; i < whole.chain.size(); i++) {
            EXPECT_EQ(windowed.chain[i].genepos, whole.chain[i].genepos) << what << " link " << i;
            EXPECT_EQ(windowed.chain[i].readpos, whole.chain[i].readpos) << what << " link " << i;
            EXPECT_EQ(windowed.chain[i].length, whole.chain[i].length) << what << " link " << i;
        }
    }

    // The read: `length` bases from `start` of the gene (start may be negative or run past the end: the
    // bases outside the gene are random), with the given positions changed.
    std::string ReadAt(std::string const& gene, long start, size_t length, std::vector<size_t> const& mismatches, std::mt19937& rng) {
        std::string read;
        for (size_t i = 0; i < length; i++) {
            long const g = start + static_cast<long>(i);
            read += (g >= 0 && g < static_cast<long>(gene.size())) ? gene[g] : "ACGT"[rng() % 4];
        }
        for (size_t m : mismatches) {
            if (m < read.size()) read[m] = read[m] == 'A' ? 'C' : 'A';
        }
        return read;
    }
}

TEST(GeneWindows, ChainsWithAnIndelBetweenTheirLinksReadBothDiagonals) {
    Reference ref({ 1945, 2400, 4500 });
    Seedmap map;
    KmerLookupSM lookup(map, 16);
    Finder finder(lookup, 15, 4, 10, *ref.loader);
    std::mt19937 rng(8);
    size_t checked = 0;
    for (uint32_t gene = 1; gene <= 3; gene++) {
        std::string const& g = ref.genes[gene - 1];
        for (long start : { 0L, 100L, 1646L, static_cast<long>(g.size()) - 153 }) {
            if (start + 153 > static_cast<long>(g.size())) continue;
            for (long deletion : { 1L, 3L, 6L }) {
                // 62 bases, a deletion of `deletion` gene bases, then the rest of the 150
                std::string read = g.substr(start, 62) + g.substr(start + 62 + deletion, 88);
                for (size_t second : { 62u, 70u, 90u }) {
                    ChainList chain = { ChainLink(static_cast<uint32_t>(start + 20), 20, 20),
                                        ChainLink(static_cast<uint32_t>(start + second + deletion), static_cast<uint16_t>(second), 20) };
                    ExpectSameExtension(finder, ref, gene, read, chain, "two links, deletion");
                    // the reverse: an insertion in the read
                    std::string inserted = g.substr(start, 62) + "ACG" + g.substr(start + 62, 85);
                    ChainList chain2 = { ChainLink(static_cast<uint32_t>(start + 20), 20, 20),
                                         ChainLink(static_cast<uint32_t>(start + second), static_cast<uint16_t>(second + 3), 20) };
                    ExpectSameExtension(finder, ref, gene, inserted, chain2, "two links, insertion");
                    checked++;
                }
            }
        }
    }
    EXPECT_GT(checked, 20u);
}

namespace {
    // Fills the undecoded part of every window with a byte for the test.
    struct FillOutside {
        explicit FillOutside(int byte) { packed::OutsideFill().store(byte); }
        ~FillOutside() { packed::OutsideFill().store(-1); }
    };
}

// One seed of a read inside the gene, over either end or (genes of 100 and 151 bases) longer than the gene, with
// mismatches or without: the finder's window extends it as the whole gene does. What lies outside a window is
// whatever the stack held; filled with a fixed byte, a read there that matters shows (and in an AddressSanitizer
// build any read there is an error). NUL is the byte that matters most to the extension loop, which stops when the
// read's NUL meets a base that is not one.
TEST(GeneWindows, NothingReadsOutsideTheWindowWhateverIsThere) {
    Reference ref({ 1500, 100, 6000, 300, 151, 1945 });
    Seedmap map;
    KmerLookupSM lookup(map, 16);
    Finder finder(lookup, 15, 4, 10, *ref.loader);
    for (int byte : std::vector<int>{ 0, '#', 'A', 'N' }) {
        FillOutside fill(byte);
        SCOPED_TRACE("outside filled with " + std::to_string(byte));
        std::mt19937 rng(9);
        size_t checked = 0;
        for (uint32_t gene = 1; gene <= ref.genes.size(); gene++) {
            std::string const& g = ref.genes[gene - 1];
            long const length = static_cast<long>(g.size());
            for (long start : { 0L, 1L, 5L, 37L, 45L, length / 3, length / 2 - 75, length - 150, length - 149, length - 100, length - 20,
                                length - 10, -20L, -30L, -140L }) {
                if (start < -140 || start >= length) continue;
                for (std::vector<size_t> mismatches : { std::vector<size_t>{}, std::vector<size_t>{ 60 }, std::vector<size_t>{ 30, 90 },
                                                        std::vector<size_t>{ 3, 40, 75, 110, 146 } }) {
                    std::string const read = ReadAt(g, start, 150, mismatches, rng);
                    // a seed of 15 bases at several places of the read, on the read's diagonal, where it lies on the gene
                    for (size_t readpos : { 0u, 20u, 40u, 62u, 100u, 120u, 135u }) {
                        long const genepos = start + static_cast<long>(readpos);
                        if (genepos < 0 || genepos + 15 > length || g.compare(genepos, 15, read, readpos, 15) != 0) continue;
                        ExpectSameExtension(finder, ref, gene, read, { ChainLink(static_cast<uint32_t>(genepos), static_cast<uint16_t>(readpos), 15) }, "one seed");
                        checked++;
                    }
                }
            }
        }
        EXPECT_GT(checked, 150u);
    }
}

namespace {
    // The maximal exact runs of at least 15 bases of the read on the diagonal (gene position - read position).
    ChainList Runs(std::string const& read, std::string const& gene, long diagonal) {
        ChainList runs;
        size_t i = 0;
        while (i < read.size()) {
            long const g = diagonal + static_cast<long>(i);
            if (g < 0 || g >= static_cast<long>(gene.size()) || read[i] != gene[g]) { i++; continue; }
            size_t j = i;
            while (j < read.size() && diagonal + static_cast<long>(j) < static_cast<long>(gene.size()) && read[j] == gene[diagonal + j]) j++;
            if (j - i >= 15) runs.emplace_back(static_cast<uint32_t>(diagonal + static_cast<long>(i)), static_cast<uint16_t>(i), static_cast<uint16_t>(j - i));
            i = j;
        }
        return runs;
    }

    // What the handler makes of random reads of random genes, as a string per anchor.
    std::vector<std::string> Alignments(Reference& ref, int fill_byte) {
        FillOutside fill(fill_byte);
        WFA2Wrapper2 aligner{4, 6, 2, 1000};
        SimpleAlignmentHandler handler(*ref.loader, aligner, 31, 3, 0.9, false);
        handler.SetAnchoredAlignment(true);
        std::mt19937 rng(10);
        std::vector<std::string> out;
        for (int round = 0; round < 400; round++) {
            uint32_t const gene = 1 + rng() % ref.genes.size();
            std::string const& g = ref.genes[gene - 1];
            long const length = static_cast<long>(g.size());
            long const start = -20 + static_cast<long>(rng() % (length + 20));  // some reads run over either end
            std::vector<size_t> mismatches;
            for (int m = 0, n = rng() % 5; m < n; m++) mismatches.push_back(rng() % 150);
            std::string read = ReadAt(g, start, 150, mismatches, rng);
            if (rng() % 6 == 0) read.erase(40 + rng() % 40, 1 + rng() % 4);  // a deletion in the read: a gap in the alignment
            ChainList chain = Runs(read, g, start);
            if (chain.empty()) continue;
            ChainAlignmentAnchor anchor(1, gene, true);
            anchor.chain = chain;
            AlignmentResult result;
            std::string rev = KmerUtils::ReverseComplement(read), id = "r";
            bool const aligned = handler.AlignAnchor(anchor, result, read, rev, false, id);
            out.push_back(std::to_string(aligned) + ' ' + (aligned ? std::to_string(result.GetAlignmentInfo().gene_alignment_start) + ' ' +
                                                                    result.GetAlignmentInfo().cigar : std::string()));
        }
        return out;
    }
}

TEST(GeneWindows, TheAlignmentOfAnAnchorDoesNotDependOnWhatIsOutsideItsWindow) {
    Reference ref({ 1500, 1945, 6000, 151, 90 });
    auto const reference = Alignments(ref, -1);
    ASSERT_GT(reference.size(), 250u);
    size_t aligned = 0;
    for (auto const& a : reference) aligned += a[0] == '1';
    EXPECT_GT(aligned, 150u);
    for (int byte : std::vector<int>{ 0, '#', 'N' }) {
        EXPECT_EQ(Alignments(ref, byte), reference) << "outside filled with " << byte;
    }
}
