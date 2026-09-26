// Unit tests for the SAM round-trip: read names written for a pair, and reading pairs, orphan
// mates and legacy files back with SamReader.
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include "IO/AlignmentOutputHandler.h"

using namespace protal;

namespace {
    // A one-gene reference (taxid 1, gene 1) written to a temporary directory and loaded.
    struct TinyReference {
        std::string gene = "ACGTTGCAAGGCTTACCGATGACTGAAACCGGTTTACGATCGGTAGCATG";  // 50 bp
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;

        TinyReference() {
            dir = std::filesystem::temp_directory_path() / ("protal_samtest_" + std::to_string(::getpid()));
            std::filesystem::create_directories(dir);
            std::string header = ">1_1\n";
            { std::ofstream fna(dir / "reference.fna"); fna << header << gene << '\n'; }
            { std::ofstream map(dir / "reference.map"); map << "1\t1\t" << header.size() << '\t' << header.size() + gene.size() << '\n'; }
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
            loader->LoadAllGenomes();
        }
        ~TinyReference() { std::filesystem::remove_all(dir); }
    };

    // An ungapped alignment of `length` bases starting at gene position `start`.
    AlignmentResult Aligned(size_t start, size_t length, bool forward) {
        AlignmentResult ar(0, 1, 1, static_cast<int32_t>(start), forward);
        auto& info = ar.GetAlignmentInfo();
        info.cigar = std::string(length, 'M');
        info.compressed_cigar = std::to_string(length) + "M";
        info.gene_alignment_start = static_cast<int>(start);
        info.alignment_length = length;
        info.matches = length;
        return ar;
    }

    FastxRecord Record(std::string id, std::string seq) {
        FastxRecord r;
        r.id = std::move(id);
        r.header = "@" + r.id;
        r.quality = std::string(seq.size(), 'I');
        r.sequence = std::move(seq);
        return r;
    }

    std::vector<std::vector<std::string>> WritePairs(TinyReference& ref, PairedAlignmentResultList results,
                                                     FastxRecord r1, FastxRecord r2) {
        std::ostringstream os;
        {
            ProtalPairedOutputHandler<false> handler(os, 5, 0, 1 << 16, *ref.loader);
            handler(results, r1, r2);
        }  // the destructor flushes the buffer
        std::vector<std::vector<std::string>> records;
        std::istringstream in(os.str());
        std::string line;
        std::vector<std::string> tokens;
        while (std::getline(in, line)) {
            LineSplitter::Split(line, "\t", tokens);
            records.push_back(tokens);
        }
        return records;
    }

    constexpr FLAG_t PAIRED = 0x1, BOTH_ALIGN = 0x2, MATE_UNMAPPED = 0x8, READ1 = 0x40, READ2 = 0x80;

    std::string SamLine(std::string const& qname, FLAG_t flag) {
        return qname + '\t' + std::to_string(flag) + "\t1_1\t1\t60\t4M\t*\t0\t4\tACGT\tIIII\tZU:i:1\tZT:i:0\n";
    }

    struct Reader {
        std::istringstream in;
        SamReader reader{ in };
        SamEntry sam1, sam2;
        bool has1 = false, has2 = false;

        explicit Reader(std::string text) : in(std::move(text)) {}
        bool Next() { return reader.Next(sam1, sam2, has1, has2); }
    };
}

TEST(PairQName, StripsOnlyAMateSuffix) {
    EXPECT_EQ(PairQName("read7/1", "read7/2"), "read7");
    EXPECT_EQ(PairQName("SRR123.1045.1", "SRR123.1045.2"), "SRR123.1045");
    EXPECT_EQ(PairQName("frag_1", "frag_2"), "frag");
    // Identical ids (Casava 1.8+, SRA without --readids) are kept whole: no neighbour collisions.
    EXPECT_EQ(PairQName("SRR123.1045", "SRR123.1045"), "SRR123.1045");
    EXPECT_EQ(PairQName("A00:8:H5:1:1101:1000:2000", "A00:8:H5:1:1101:1000:2000"), "A00:8:H5:1:1101:1000:2000");
    // Unrelated ids fall back to the first mate's id for both mates.
    EXPECT_EQ(PairQName("readA", "readB"), "readA");
}

TEST(ReadQName, StripsSlashMateSuffixOnly) {
    EXPECT_EQ(ReadQName("read7/1"), "read7");
    EXPECT_EQ(ReadQName("read7/2"), "read7");
    EXPECT_EQ(ReadQName("SRR123.1045"), "SRR123.1045");
    EXPECT_EQ(ReadQName("x"), "x");
}

TEST(SamReader, ReadsMatesTogether) {
    Reader r(SamLine("a", PAIRED | BOTH_ALIGN | READ1) + SamLine("a", PAIRED | BOTH_ALIGN | READ2));
    ASSERT_TRUE(r.Next());
    EXPECT_TRUE(r.has1);
    EXPECT_TRUE(r.has2);
    EXPECT_EQ(r.sam1.m_qname, "a");
    EXPECT_EQ(r.sam2.m_qname, "a");
    EXPECT_FALSE(r.Next());
}

TEST(SamReader, FlaggedOrphanRead1DoesNotSwallowNextRead) {
    Reader r(SamLine("orphan", PAIRED | READ1 | MATE_UNMAPPED) +
             SamLine("b", PAIRED | BOTH_ALIGN | READ1) + SamLine("b", PAIRED | BOTH_ALIGN | READ2));
    ASSERT_TRUE(r.Next());
    EXPECT_TRUE(r.has1);
    EXPECT_FALSE(r.has2);
    EXPECT_EQ(r.sam1.m_qname, "orphan");

    ASSERT_TRUE(r.Next());
    EXPECT_TRUE(r.has1 && r.has2);
    EXPECT_EQ(r.sam1.m_qname, "b");
    EXPECT_EQ(r.sam2.m_qname, "b");
    EXPECT_FALSE(r.Next());
}

TEST(SamReader, LegacyUnflaggedOrphanKeepsTheLookAheadRecord) {
    // Older protal versions wrote an orphan read1 with 0x1 but without 0x8.
    Reader r(SamLine("legacy", PAIRED | READ1) +
             SamLine("c", PAIRED | BOTH_ALIGN | READ1) + SamLine("c", PAIRED | BOTH_ALIGN | READ2));
    ASSERT_TRUE(r.Next());
    EXPECT_TRUE(r.has1);
    EXPECT_FALSE(r.has2);
    EXPECT_EQ(r.sam1.m_qname, "legacy");

    ASSERT_TRUE(r.Next());
    EXPECT_TRUE(r.has1 && r.has2);
    EXPECT_EQ(r.sam1.m_qname, "c");
    EXPECT_FALSE(r.Next());
}

TEST(SamReader, Read2OnlyAndTrailingOrphan) {
    Reader r(SamLine("d", PAIRED | READ2 | MATE_UNMAPPED) + SamLine("e", PAIRED | READ1));
    ASSERT_TRUE(r.Next());
    EXPECT_FALSE(r.has1);
    EXPECT_TRUE(r.has2);
    EXPECT_EQ(r.sam2.m_qname, "d");

    // A paired read1 on the last line has no mate to read: it is still returned.
    ASSERT_TRUE(r.Next());
    EXPECT_TRUE(r.has1);
    EXPECT_FALSE(r.has2);
    EXPECT_EQ(r.sam1.m_qname, "e");
    EXPECT_FALSE(r.Next());
}

TEST(PairedOutputHandler, WritesPairsWhereOnlyRead1Aligned) {
    TinyReference ref;
    auto r1 = Record("frag/1", ref.gene.substr(0, 20));
    auto r2 = Record("frag/2", "NNNNNNNNNNNNNNNNNNNN");
    auto records = WritePairs(ref, { { Aligned(0, 20, true), AlignmentResult() } }, r1, r2);

    ASSERT_EQ(records.size(), 1u);
    FLAG_t flag = std::stoul(records[0][1]);
    EXPECT_EQ(records[0][0], "frag");
    EXPECT_TRUE(Flag::IsRead1(flag));
    EXPECT_TRUE(Flag::IsMateUnmapped(flag));
    EXPECT_FALSE(Flag::IsPairBothAlign(flag));
}

TEST(PairedOutputHandler, WritesBothMatesAdjacentWithOneName) {
    TinyReference ref;
    auto r1 = Record("frag/1", ref.gene.substr(0, 20));
    auto r2 = Record("frag/2", KmerUtils::ReverseComplement(ref.gene.substr(24, 20)));
    auto records = WritePairs(ref, { { Aligned(0, 20, true), Aligned(24, 20, false) } }, r1, r2);

    ASSERT_EQ(records.size(), 2u);
    EXPECT_EQ(records[0][0], "frag");
    EXPECT_EQ(records[1][0], "frag");
    FLAG_t flag1 = std::stoul(records[0][1]);
    FLAG_t flag2 = std::stoul(records[1][1]);
    EXPECT_TRUE(Flag::IsRead1(flag1) && Flag::IsPairBothAlign(flag1) && !Flag::IsMateUnmapped(flag1));
    EXPECT_TRUE(Flag::IsRead2(flag2) && Flag::IsReverseComplement(flag2));
    EXPECT_EQ(records[1][9], ref.gene.substr(24, 20));  // SEQ is stored in reference orientation
}

TEST(PairedOutputHandler, SkipsOnlyTheInconsistentCandidate) {
    TinyReference ref;
    auto r1 = Record("frag/1", ref.gene.substr(0, 20));
    auto r2 = Record("frag/2", "NNNNNNNNNNNNNNNNNNNN");
    // The first candidate places the read at the wrong position, so its CIGAR disagrees with the
    // sequences; the second candidate is correct and must still be written, as the primary.
    auto records = WritePairs(ref, { { Aligned(3, 20, true), AlignmentResult() },
                                     { Aligned(0, 20, true), AlignmentResult() } }, r1, r2);

    ASSERT_EQ(records.size(), 1u);
    EXPECT_EQ(records[0][3], "1");  // POS of the correct candidate (1-based)
    EXPECT_FALSE(Flag::IsNotPrimaryAlignment(std::stoul(records[0][1])));
}
