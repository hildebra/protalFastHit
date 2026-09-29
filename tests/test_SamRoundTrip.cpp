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
    // A two-gene reference (taxid 1, genes 1 and 2) written to a temporary directory and loaded.
    struct TinyReference {
        std::string gene = "ACGTTGCAAGGCTTACCGATGACTGAAACCGGTTTACGATCGGTAGCATG";   // 50 bp, gene 1
        std::string gene2 = "TTGACCAGTCAGGATCCATTGCAGGTACTTGACCGTAAGCTGCATTGACA";  // 50 bp, gene 2
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;

        TinyReference() {
            dir = std::filesystem::temp_directory_path() / ("protal_samtest_" + std::to_string(::getpid()));
            std::filesystem::create_directories(dir);
            std::ofstream fna(dir / "reference.fna"), map(dir / "reference.map");
            size_t offset = 0;
            for (auto const& [id, seq] : { std::pair{ 1, gene }, std::pair{ 2, gene2 } }) {
                std::string header = ">1_" + std::to_string(id) + "\n";
                fna << header << seq << '\n';
                map << "1\t" << id << '\t' << offset + header.size() << '\t' << offset + header.size() + seq.size() << '\n';
                offset += header.size() + seq.size() + 1;
            }
            fna.close();
            map.close();
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
            loader->LoadAllGenomes();
        }
        ~TinyReference() { std::filesystem::remove_all(dir); }
    };

    // An ungapped alignment of `length` bases starting at position `start` of gene `geneid`.
    AlignmentResult Aligned(size_t start, size_t length, bool forward, uint32_t geneid = 1) {
        AlignmentResult ar(0, 1, geneid, static_cast<int32_t>(start), forward);
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

TEST(PairedOutputHandler, WritesPairsWhereOnlyRead2Aligned) {
    TinyReference ref;
    auto r1 = Record("frag/1", "NNNNNNNNNNNNNNNNNNNN");
    auto r2 = Record("frag/2", KmerUtils::ReverseComplement(ref.gene.substr(10, 20)));
    auto records = WritePairs(ref, { { AlignmentResult(), Aligned(10, 20, false) } }, r1, r2);

    ASSERT_EQ(records.size(), 1u);
    FLAG_t flag = std::stoul(records[0][1]);
    EXPECT_TRUE(Flag::IsRead2(flag));
    EXPECT_TRUE(Flag::IsMateUnmapped(flag));
    EXPECT_EQ(records[0][3], "11");
    EXPECT_GT(std::stoi(records[0][4]), 4);
}

TEST(PairedOutputHandler, WritesMatesOnTwoGenesWithTheirOwnMapq) {
    // The fragment spans the end of gene 1 and the start of gene 2: candidates pair mates only
    // within one gene, so each mate is a candidate of its own.
    TinyReference ref;
    auto r1 = Record("frag/1", ref.gene.substr(30, 20));
    auto r2 = Record("frag/2", KmerUtils::ReverseComplement(ref.gene2.substr(0, 20)));
    auto records = WritePairs(ref, { { Aligned(30, 20, true), AlignmentResult() },
                                     { AlignmentResult(), Aligned(0, 20, false, 2) } }, r1, r2);

    ASSERT_EQ(records.size(), 2u);
    FLAG_t flag1 = std::stoul(records[0][1]);
    FLAG_t flag2 = std::stoul(records[1][1]);
    EXPECT_TRUE(Flag::IsRead1(flag1) && Flag::IsRead2(flag2));
    for (FLAG_t flag : { flag1, flag2 }) {
        EXPECT_TRUE(Flag::IsPaired(flag));
        EXPECT_FALSE(Flag::IsPairBothAlign(flag));
        EXPECT_FALSE(Flag::IsMateUnmapped(flag));
        EXPECT_FALSE(Flag::IsNotPrimaryAlignment(flag));
    }
    EXPECT_EQ(records[0][2], "1_1");
    EXPECT_EQ(records[1][2], "1_2");
    EXPECT_EQ(records[0][6], "1_2");  // RNEXT
    EXPECT_EQ(records[1][6], "1_1");
    EXPECT_GT(std::stoi(records[0][4]), 4);  // each mate is unambiguous on its own
    EXPECT_GT(std::stoi(records[1][4]), 4);

    // Read back, they are one fragment with a mate on each gene.
    std::string text;
    for (auto const& record : records) {
        for (size_t i = 0; i < record.size(); i++) text += (i ? "\t" : "") + record[i];
        text += '\n';
    }
    std::istringstream in(text);
    SamReader reader(in);
    SamEntry sam1, sam2;
    bool has1 = false, has2 = false;
    ASSERT_TRUE(reader.Next(sam1, sam2, has1, has2));
    EXPECT_TRUE(has1 && has2);
    EXPECT_FALSE(reader.Next(sam1, sam2, has1, has2));
}

TEST(PairedOutputHandler, SplitMatesOnOneGeneKeepTheirPositions) {
    // Both mates on gene 1 but in an orientation that does not form a pair.
    TinyReference ref;
    auto r1 = Record("frag/1", ref.gene.substr(0, 20));
    auto r2 = Record("frag/2", ref.gene.substr(25, 20));
    auto records = WritePairs(ref, { { Aligned(0, 20, true), AlignmentResult() },
                                     { AlignmentResult(), Aligned(25, 20, true) } }, r1, r2);
    ASSERT_EQ(records.size(), 2u);
    EXPECT_EQ(records[0][6], "=");
    EXPECT_EQ(records[0][7], "26");  // PNEXT
    EXPECT_EQ(records[1][7], "1");
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

TEST(SamReader, SingleEndReadsAreFirstReads) {
    Reader r(SamLine("fwd", 0) + SamLine("rev", 0x10) + SamLine("rev", 0x10 | 0x100));
    for (std::string const name : { "fwd", "rev", "rev" }) {
        ASSERT_TRUE(r.Next());
        EXPECT_TRUE(r.has1);
        EXPECT_FALSE(r.has2);
        EXPECT_EQ(r.sam1.m_qname, name);
    }
    EXPECT_FALSE(r.Next());
    EXPECT_EQ(r.reader.Records(), 3u);
    EXPECT_EQ(r.reader.PairedRecords(), 0u);
}

TEST(SamReader, TellsPairedFromSingleEndReads) {
    auto paired = [](std::string text) {
        std::istringstream in(std::move(text));
        return HoldsPairedReads(in);
    };
    EXPECT_EQ(paired("@HD\tVN:1.6\n" + SamLine("a", PAIRED | BOTH_ALIGN | READ1) + SamLine("a", PAIRED | BOTH_ALIGN | READ2)), true);
    EXPECT_EQ(paired(SamLine("orphan", PAIRED | READ2 | MATE_UNMAPPED)), true);
    EXPECT_EQ(paired("@HD\tVN:1.6\n" + SamLine("a", 0x10) + SamLine("b", 0)), false);
    EXPECT_EQ(paired("@HD\tVN:1.6\n@SQ\tSN:1_1\tLN:50\n"), std::nullopt);
}

namespace {
    std::vector<std::vector<std::string>> WriteSingle(TinyReference& ref, AlignmentResultList results, FastxRecord record,
                                                      size_t max_out = 5) {
        std::ostringstream os;
        {
            ProtalSingleOutputHandler<false> handler(os, max_out, 0, 1 << 16, *ref.loader);
            handler(results, record);
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

    // An alignment of `length` bases with `mismatches` of them counted as mismatches (for the score
    // only; the CIGAR stays all M).
    AlignmentResult Weaker(size_t start, size_t length, size_t mismatches, uint32_t geneid) {
        auto ar = Aligned(start, length, true, geneid);
        ar.GetAlignmentInfo().matches = length - mismatches;
        ar.GetAlignmentInfo().mismatches = mismatches;
        return ar;
    }
}

TEST(SingleOutputHandler, WritesUnpairedRecordsInReferenceOrientation) {
    TinyReference ref;
    auto fwd = WriteSingle(ref, { Aligned(5, 20, true) }, Record("read1/1", ref.gene.substr(5, 20)));
    ASSERT_EQ(fwd.size(), 1u);
    EXPECT_EQ(fwd[0][0], "read1");  // QNAME without a /1 suffix
    EXPECT_EQ(std::stoul(fwd[0][1]), 0u);  // not paired, forward, primary
    EXPECT_EQ(fwd[0][2], "1_1");
    EXPECT_EQ(fwd[0][3], "6");
    EXPECT_GT(std::stoi(fwd[0][4]), 4);  // the only candidate: unambiguous
    EXPECT_EQ(fwd[0][6], "*");
    EXPECT_EQ(fwd[0][7], "0");
    EXPECT_EQ(fwd[0][8], "0");  // TLEN

    auto rev = WriteSingle(ref, { Aligned(10, 20, false) }, Record("read2", KmerUtils::ReverseComplement(ref.gene.substr(10, 20))));
    ASSERT_EQ(rev.size(), 1u);
    FLAG_t flag = std::stoul(rev[0][1]);
    EXPECT_TRUE(Flag::IsReverseComplement(flag));
    EXPECT_FALSE(Flag::IsPaired(flag) || Flag::IsRead1(flag) || Flag::IsRead2(flag) || Flag::IsUnmapped(flag));
    EXPECT_EQ(rev[0][9], ref.gene.substr(10, 20));  // SEQ is stored in reference orientation

    // Read back, each record is a single read in the first slot.
    std::string text;
    for (auto const& record : { fwd[0], rev[0] }) {
        for (size_t i = 0; i < record.size(); i++) text += (i ? "\t" : "") + record[i];
        text += '\n';
    }
    Reader r(text);
    ASSERT_TRUE(r.Next());
    EXPECT_TRUE(r.has1 && !r.has2);
    EXPECT_EQ(r.sam1.m_qname, "read1");
    ASSERT_TRUE(r.Next());
    EXPECT_TRUE(r.has1 && !r.has2);
    EXPECT_TRUE(r.sam1.IsReversed());
    EXPECT_FALSE(r.Next());
}

TEST(SingleOutputHandler, RanksCandidatesAndWritesTheOthersAsSecondary) {
    TinyReference ref;
    // A read of Ns matches any gene position, so both candidates are consistent.
    auto records = WriteSingle(ref, { Weaker(3, 20, 5, 2), Aligned(0, 20, true) }, Record("r", std::string(20, 'N')));
    ASSERT_EQ(records.size(), 2u);
    EXPECT_EQ(records[0][2], "1_1");  // the better candidate first
    EXPECT_FALSE(Flag::IsNotPrimaryAlignment(std::stoul(records[0][1])));
    EXPECT_EQ(std::stoi(records[0][4]), MAPQv2(40, 15));  // bitscores 20*2 and 15*2-5*3
    EXPECT_EQ(records[1][2], "1_2");
    EXPECT_TRUE(Flag::IsNotPrimaryAlignment(std::stoul(records[1][1])));
    EXPECT_EQ(records[1][4], "0");
    EXPECT_EQ(records[0][0], records[1][0]);

    auto first_only = WriteSingle(ref, { Weaker(3, 20, 5, 2), Aligned(0, 20, true) }, Record("r", std::string(20, 'N')), 1);
    ASSERT_EQ(first_only.size(), 1u);
    EXPECT_EQ(first_only[0][2], "1_1");
}

TEST(SingleOutputHandler, TheSameAlignmentTwiceIsOneCandidate) {
    // Two anchors of one read on one gene can give the same alignment; that is no second hit.
    TinyReference ref;
    auto records = WriteSingle(ref, { Aligned(0, 20, true), Aligned(0, 20, true) }, Record("r", ref.gene.substr(0, 20)));
    ASSERT_EQ(records.size(), 1u);
    EXPECT_EQ(std::stoi(records[0][4]), MAPQv2(40, 0));
}

TEST(SingleOutputHandler, SkipsOnlyTheInconsistentCandidate) {
    TinyReference ref;
    auto records = WriteSingle(ref, { Aligned(3, 20, true), Aligned(0, 20, true) }, Record("r", ref.gene.substr(0, 20)));
    ASSERT_EQ(records.size(), 1u);
    EXPECT_EQ(records[0][3], "1");
    EXPECT_FALSE(Flag::IsNotPrimaryAlignment(std::stoul(records[0][1])));
}
