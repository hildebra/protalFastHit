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
#include "Core/MateGuidance.h"

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

    // The records written; `genes`, if given, receives the genes the handler reported for them.
    std::vector<std::vector<std::string>> WritePairs(TinyReference& ref, PairedAlignmentResultList results,
                                                     FastxRecord r1, FastxRecord r2, std::vector<uint64_t>* genes = nullptr) {
        std::ostringstream os;
        SamStreamSink sink(os);
        {
            ProtalPairedOutputHandler<false> handler(sink, 5, 0, 1 << 16, *ref.loader);
            handler(results, r1, r2);
        }  // the destructor flushes the buffer
        if (genes) *genes = sink.Genes();
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
    std::vector<uint64_t> genes;
    auto records = WritePairs(ref, { { Aligned(30, 20, true), AlignmentResult() },
                                     { AlignmentResult(), Aligned(0, 20, false, 2) } }, r1, r2, &genes);
    // Both genes go into the SAM header: each is a record's RNAME and the other's RNEXT.
    EXPECT_EQ(genes, (std::vector<uint64_t>{ SamGeneKey(1, 1), SamGeneKey(1, 2) }));

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

TEST(PairedOutputHandler, ReportsOnlyTheGenesOfRecordsWritten) {
    TinyReference ref;
    auto r1 = Record("frag/1", ref.gene.substr(0, 20));
    auto r2 = Record("frag/2", "NNNNNNNNNNNNNNNNNNNN");
    // The candidate on gene 2 does not fit gene 2 and is skipped: gene 2 must not be in the header.
    std::vector<uint64_t> genes;
    auto records = WritePairs(ref, { { Aligned(0, 20, true, 2), AlignmentResult() },
                                     { Aligned(0, 20, true), AlignmentResult() } }, r1, r2, &genes);
    ASSERT_EQ(records.size(), 1u);
    EXPECT_EQ(records[0][2], "1_1");
    EXPECT_EQ(genes, (std::vector<uint64_t>{ SamGeneKey(1, 1) }));
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
                                                      size_t max_out = 5, std::vector<uint64_t>* genes = nullptr) {
        std::ostringstream os;
        SamStreamSink sink(os);
        {
            ProtalSingleOutputHandler<false> handler(sink, max_out, 0, 1 << 16, *ref.loader);
            handler(results, record);
        }  // the destructor flushes the buffer
        if (genes) *genes = sink.Genes();
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

    std::vector<uint64_t> genes;
    auto first_only = WriteSingle(ref, { Weaker(3, 20, 5, 2), Aligned(0, 20, true) }, Record("r", std::string(20, 'N')), 1, &genes);
    ASSERT_EQ(first_only.size(), 1u);
    EXPECT_EQ(first_only[0][2], "1_1");
    EXPECT_EQ(genes, (std::vector<uint64_t>{ SamGeneKey(1, 1) }));  // gene 2's candidate was not written
    WriteSingle(ref, { Weaker(3, 20, 5, 2), Aligned(0, 20, true) }, Record("r", std::string(20, 'N')), 5, &genes);
    EXPECT_EQ(genes, (std::vector<uint64_t>{ SamGeneKey(1, 1), SamGeneKey(1, 2) }));
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

// A read that seeded on taxa but aligned nowhere: counted per taxon for the SAM header by default, an unmapped record
// with --write_unmapped_reads; the reader counts the header line as it counts the records.
TEST(SingleOutputHandler, CountsUnalignedReadsForTheHeaderOrWritesTheirRecords) {
    TinyReference ref;
    auto failed_of = [](SamReader& reader) {
        SamEntry a, b;
        bool has_a = false, has_b = false;
        EXPECT_FALSE(reader.Next(a, b, has_a, has_b));
        auto counts = reader.FailedCandidates();
        counts.resize(8, 0);
        return counts;
    };
    std::vector<uint32_t> const expected{ 0, 0, 0, 1, 0, 0, 0, 2 };
    for (bool const write_records : { false, true }) {
        std::ostringstream os;
        SamStreamSink sink(os);
        sink.SetUnmappedRecords(write_records);
        {
            ProtalSingleOutputHandler<false> handler(sink, 5, 0, 1 << 16, *ref.loader);
            AlignmentResultList none, none2;
            auto record = Record("read1", ref.gene.substr(5, 20));
            auto record2 = Record("read2", ref.gene.substr(5, 20));
            handler(none, record, { 3, 7 });
            handler(none2, record2, { 7 });
        }  // the destructor hands over the counts
        auto const counts = sink.FailedCandidates();
        if (write_records) {
            EXPECT_TRUE(std::all_of(counts.begin(), counts.end(), [](uint64_t n) { return n == 0; }));
            std::istringstream in(os.str());
            SamReader reader(in);
            EXPECT_EQ(failed_of(reader), expected);
            EXPECT_EQ(reader.Skipped().at("unmapped"), 2u);
        } else {
            EXPECT_EQ(os.str(), "");
            EXPECT_EQ(FailedCandidatesLine(counts), kSamFailedCandidatesComment + "3:1,7:2\n");
            std::istringstream in("@HD\tVN:1.6\n" + FailedCandidatesLine(counts));
            SamReader reader(in);
            EXPECT_EQ(failed_of(reader), expected);
            EXPECT_TRUE(reader.Skipped().empty());
        }
    }
}

TEST(AlternativesTag, ListsOtherTaxaByTheirBestCandidate) {
    EXPECT_EQ(AlignmentEdits("2S19M1X128M"), 3);  // clipped bases count as edits
    EXPECT_EQ(AlignmentEdits("5M2I3M1D4M"), 3);
    EXPECT_EQ(AlignmentEdits("150M"), 0);
    // The best is taxon 1 with 2 edits: its own other candidates are not listed, each other taxon once with its
    // best candidate, fewest edits first, none more than kAlternativeMaxEdits worse.
    EXPECT_EQ(AlternativesTag(1, 2, { { 1, 2 }, { 1, 0 }, { 2, 5 }, { 2, 3 }, { 3, 2 }, { 4, 9 } }), "3:0,2:1");
    EXPECT_EQ(AlternativesTag(1, 2, { { 1, 2 } }), "*");
    EXPECT_EQ(AlternativesTag(1, 4, { { 2, 3 } }), "2:-1");  // fewer edits but a lower score: still an alternative
    std::vector<std::pair<uint32_t, int>> many;
    for (uint32_t t = 2; t < 10; t++) many.emplace_back(t, 2);
    EXPECT_EQ(AlternativesTag(1, 2, many), "2:0,3:0,4:0,5:0");  // at most kAlternativesListed
}

namespace {
    // A candidate of another taxon at the start of gene 1, `edits` of its 20 bases mismatches (in the CIGAR and
    // the score).
    AlignmentResult OfTaxon(uint32_t taxid, size_t edits) {
        AlignmentResult ar(0, taxid, 1, 0, true);
        auto& info = ar.GetAlignmentInfo();
        info.cigar = std::string(20 - edits, 'M') + std::string(edits, 'X');
        info.compressed_cigar = std::to_string(20 - edits) + "M" + (edits ? std::to_string(edits) + "X" : "");
        info.gene_alignment_start = 0;
        info.alignment_length = 20;
        info.matches = 20 - edits;
        info.mismatches = edits;
        return ar;
    }

    std::string TagOf(std::vector<std::string> const& record, std::string const& name) {
        for (size_t i = 11; i < record.size(); i++) {
            if (record[i].compare(0, name.size() + 3, name + ":Z:") == 0) return record[i].substr(name.size() + 3);
        }
        return "absent";
    }
}

TEST(SingleOutputHandler, TagsTheBestRecordWithTheReadsAlternativesInOtherTaxa) {
    TinyReference ref;
    // With -m 1 only the best is written (taxon 2's candidate is not, so its genome is not needed).
    auto records = WriteSingle(ref, { OfTaxon(2, 1), Aligned(0, 20, true), Weaker(3, 20, 5, 2) },
                               Record("r", ref.gene.substr(0, 20)), 1);
    ASSERT_EQ(records.size(), 1u);
    EXPECT_EQ(records[0][2], "1_1");
    EXPECT_EQ(TagOf(records[0], "ZA"), "2:1");  // taxon 1's own other gene is no alternative

    auto alone = WriteSingle(ref, { Aligned(0, 20, true) }, Record("r", ref.gene.substr(0, 20)));
    ASSERT_EQ(alone.size(), 1u);
    EXPECT_EQ(TagOf(alone[0], "ZA"), "*");

    // Secondary records carry no ZA.
    auto both = WriteSingle(ref, { Weaker(3, 20, 5, 2), Aligned(0, 20, true) }, Record("r", std::string(20, 'N')));
    ASSERT_EQ(both.size(), 2u);
    EXPECT_EQ(TagOf(both[0], "ZA"), "*");
    EXPECT_EQ(TagOf(both[1], "ZA"), "absent");
}

namespace {
    // Taxa 1 and 2, each with genes 1 and 2 of TinyReference's sequences.
    struct TwoTaxaReference {
        TinyReference tiny;
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;

        TwoTaxaReference() {
            dir = tiny.dir / "two";
            std::filesystem::create_directories(dir);
            std::ofstream fna(dir / "reference.fna"), map(dir / "reference.map");
            size_t offset = 0;
            for (int taxid : { 1, 2 }) {
                for (auto const& [id, seq] : { std::pair{ 1, tiny.gene }, std::pair{ 2, tiny.gene2 } }) {
                    std::string header = ">" + std::to_string(taxid) + "_" + std::to_string(id) + "\n";
                    fna << header << seq << '\n';
                    map << taxid << '\t' << id << '\t' << offset + header.size() << '\t' << offset + header.size() + seq.size() << '\n';
                    offset += header.size() + seq.size() + 1;
                }
            }
            fna.close();
            map.close();
            loader = std::make_unique<GenomeLoader>((dir / "reference.fna").string(), (dir / "reference.map").string());
            loader->LoadAllGenomes();
        }
    };

    // An ungapped alignment of taxon `taxid`; `mismatches` count for the score only.
    AlignmentResult OnTaxon(uint32_t taxid, uint32_t geneid, size_t start, bool forward, size_t mismatches = 0) {
        AlignmentResult ar(0, taxid, geneid, static_cast<int32_t>(start), forward);
        auto& info = ar.GetAlignmentInfo();
        info.cigar = std::string(20, 'M');
        info.compressed_cigar = "20M";
        info.gene_alignment_start = static_cast<int>(start);
        info.alignment_length = 20;
        info.matches = 20 - mismatches;
        info.mismatches = mismatches;
        return ar;
    }
}

TEST(PairedOutputHandler, MatesOnTwoGenesTakeThePairsConsensusTaxon) {
    TwoTaxaReference ref;
    auto r1 = Record("frag/1", ref.tiny.gene.substr(30, 20));
    auto r2 = Record("frag/2", KmerUtils::ReverseComplement(ref.tiny.gene2.substr(0, 20)));
    auto write = [&](PairedAlignmentResultList results) {
        std::ostringstream os;
        SamStreamSink sink(os);
        {
            ProtalPairedOutputHandler<false> handler(sink, 5, 0, 1 << 16, *ref.loader);
            handler(results, r1, r2);
        }
        std::vector<std::vector<std::string>> records;
        std::istringstream in(os.str());
        std::string line;
        std::vector<std::string> tokens;
        while (std::getline(in, line)) {
            LineSplitter::Split(line, "\t", tokens);
            records.push_back(tokens);
        }
        return records;
    };
    // Mate 1 fits taxon 1 clearly better than taxon 2; mate 2 fits taxon 2 a little better than taxon 1. On both
    // mates together taxon 1 wins (40 + 35 against 25 + 40): mate 2 takes its alignment of taxon 1.
    auto both = write({ { OnTaxon(1, 1, 30, true), AlignmentResult() }, { OnTaxon(2, 1, 30, true, 3), AlignmentResult() },
                        { AlignmentResult(), OnTaxon(2, 2, 0, false) }, { AlignmentResult(), OnTaxon(1, 2, 0, false, 1) } });
    ASSERT_EQ(both.size(), 2u);
    EXPECT_EQ(both[0][2], "1_1");
    EXPECT_EQ(both[1][2], "1_2");
    EXPECT_EQ(std::stoi(both[1][4]), MAPQv2(75, 65));
    EXPECT_GT(std::stoi(both[0][4]), 4);

    // Mate 2 has no alignment of taxon 1 (the database lacks that gene of the species): it keeps its alignment of
    // taxon 2, but with MAPQ 0, so that it is no evidence for taxon 2.
    auto lacking = write({ { OnTaxon(1, 1, 30, true), AlignmentResult() }, { OnTaxon(2, 1, 30, true, 3), AlignmentResult() },
                           { AlignmentResult(), OnTaxon(2, 2, 0, false) } });
    ASSERT_EQ(lacking.size(), 2u);
    EXPECT_EQ(lacking[0][2], "1_1");
    EXPECT_EQ(lacking[1][2], "2_2");
    EXPECT_EQ(lacking[1][4], "0");
}

namespace {
    // Stands in for the alignment handler: records the anchors it is given and aligns each where it says, ungapped.
    struct RecordingHandler {
        std::vector<CAlignmentAnchor> seen;
        bool fail = false;
        void operator()(AlignmentAnchorList& anchors, AlignmentResultList& results, std::string const& sequence,
                        std::string const&, size_t, std::string&) {
            for (auto const& a : anchors) {
                seen.push_back(a);
                if (fail) continue;
                auto const start = static_cast<int32_t>(a.chain.front().genepos) - static_cast<int32_t>(a.chain.front().readpos);
                AlignmentResult ar(0, a.taxid, a.geneid, start, a.forward);
                auto& info = ar.GetAlignmentInfo();
                info.cigar = std::string(sequence.size(), 'M');
                info.compressed_cigar = std::to_string(sequence.size()) + "M";
                info.gene_alignment_start = start;
                info.alignment_length = sequence.size();
                info.matches = sequence.size();
                results.push_back(std::move(ar));
            }
        }
    };
}

TEST(MateGuidance, TheDiagonalOfAMatesKmers) {
    std::string const gene = "ACGTTGCAAGGCTTACCGATGACTGAAACCGGTTTACGATCGGTAGCATG";
    auto const d = mate_guidance::BestDiagonal(gene.substr(24, 20), gene, 0, 12);
    EXPECT_EQ(d.offset, 24);
    EXPECT_EQ(d.hits, 9u);  // 20 - 12 + 1 k-mers
    EXPECT_EQ(d.readpos, 0u);
    auto const shifted = mate_guidance::BestDiagonal(gene.substr(24, 20), std::string_view(gene).substr(10), 10, 12);
    EXPECT_EQ(shifted.offset, 24);  // window positions are gene positions
    EXPECT_EQ(mate_guidance::BestDiagonal(std::string(20, 'N'), gene, 0, 12).hits, 0u);
}

TEST(MateGuidance, ASureMateGuidesTheOther) {
    TinyReference ref;
    auto r1 = Record("frag/1", ref.gene.substr(0, 20));
    auto r2 = Record("frag/2", KmerUtils::ReverseComplement(ref.gene.substr(24, 20)));
    std::string const rev1 = KmerUtils::ReverseComplement(r1.sequence), rev2 = KmerUtils::ReverseComplement(r2.sequence);
    MateGuidanceCounts counts;

    // Mate 2 has no anchor: it is looked for downstream of mate 1 on gene 1 and found at 24, reverse; the best
    // candidate becomes a pair.
    {
        PairedAlignmentResultList pairs{ { Aligned(0, 20, true), AlignmentResult() } };
        RecordingHandler handler;
        EXPECT_TRUE(GuideMate(pairs, {}, {}, r1, r2, rev1, rev2, handler, *ref.loader, counts));
        ASSERT_EQ(handler.seen.size(), 1u);
        EXPECT_EQ(handler.seen[0].taxid, 1u);
        EXPECT_EQ(handler.seen[0].geneid, 1u);
        EXPECT_FALSE(handler.seen[0].forward);
        EXPECT_EQ(static_cast<int>(handler.seen[0].chain[0].genepos) - handler.seen[0].chain[0].readpos, 24);
        ASSERT_EQ(pairs.size(), 1u);
        ASSERT_TRUE(pairs[0].second.IsSet());
        EXPECT_EQ(pairs[0].second.GetAlignmentInfo().gene_alignment_start, 24);
        EXPECT_EQ(counts.rescued, 1u);
    }
    // Mate 2 has an anchor of the taxon on gene 2: that one is aligned, as a candidate of its own.
    {
        PairedAlignmentResultList pairs{ { Aligned(0, 20, true), AlignmentResult() } };
        CAlignmentAnchor anchor(1, 2, false);
        anchor.chain.emplace_back(ChainLink(5, 0, 20));
        anchor.total_length = 20;
        RecordingHandler handler;
        EXPECT_TRUE(GuideMate(pairs, {}, { anchor }, r1, r2, rev1, rev2, handler, *ref.loader, counts));
        ASSERT_EQ(handler.seen.size(), 1u);
        EXPECT_EQ(handler.seen[0].geneid, 2u);
        ASSERT_EQ(pairs.size(), 2u);
        EXPECT_FALSE(pairs[1].first.IsSet());
        EXPECT_EQ(pairs[1].second.GeneId(), 2u);
        EXPECT_EQ(counts.from_anchor, 1u);
    }
    // No guidance: the best candidate is a pair; mate 1 is not sure (a second alignment nearly as good); mate 2 has a
    // candidate of mate 1's taxon already (the consensus decides between them).
    {
        RecordingHandler handler;
        PairedAlignmentResultList paired{ { Aligned(0, 20, true), Aligned(24, 20, false) } };
        EXPECT_FALSE(GuideMate(paired, {}, {}, r1, r2, rev1, rev2, handler, *ref.loader, counts));
        PairedAlignmentResultList unsure{ { Aligned(0, 20, true), AlignmentResult() }, { Aligned(10, 20, true, 2), AlignmentResult() } };
        EXPECT_FALSE(GuideMate(unsure, {}, {}, r1, r2, rev1, rev2, handler, *ref.loader, counts));
        PairedAlignmentResultList has_taxon{ { Aligned(0, 20, true), AlignmentResult() }, { AlignmentResult(), Aligned(5, 20, false, 2) } };
        EXPECT_FALSE(GuideMate(has_taxon, {}, {}, r1, r2, rev1, rev2, handler, *ref.loader, counts));
        EXPECT_TRUE(handler.seen.empty());
    }
    // A forward mate at the gene's end: fewer than kRescueMinBases of the gene lie downstream, nothing is looked for.
    {
        PairedAlignmentResultList pairs{ { Aligned(40, 10, true), AlignmentResult() } };
        RecordingHandler handler;
        size_t const looked_for = counts.looked_for;
        EXPECT_FALSE(GuideMate(pairs, {}, {}, r1, r2, rev1, rev2, handler, *ref.loader, counts));
        EXPECT_EQ(counts.looked_for, looked_for);
        EXPECT_TRUE(handler.seen.empty());
    }
    // A reverse mate: the other is looked for upstream, forward; mate 2 does not fit the gene that way, so no k-mers
    // place it and nothing is aligned.
    {
        PairedAlignmentResultList pairs{ { Aligned(30, 20, false), AlignmentResult() } };
        RecordingHandler handler;
        size_t const looked_for = counts.looked_for;
        EXPECT_FALSE(GuideMate(pairs, {}, {}, r1, r2, rev1, rev2, handler, *ref.loader, counts));
        EXPECT_EQ(counts.looked_for, looked_for + 1);
        EXPECT_TRUE(handler.seen.empty());
    }
}

TEST(PairedOutputHandler, TagsBothMatesOfTheBestPair) {
    TinyReference ref;
    auto r1 = Record("frag/1", ref.gene.substr(0, 20));
    auto r2 = Record("frag/2", KmerUtils::ReverseComplement(ref.gene.substr(24, 20)));
    auto records = WritePairs(ref, { { Aligned(0, 20, true), Aligned(24, 20, false) } }, r1, r2);
    ASSERT_EQ(records.size(), 2u);
    EXPECT_EQ(TagOf(records[0], "ZA"), "*");
    EXPECT_EQ(TagOf(records[1], "ZA"), "*");

    // Mates on two genes (written as split mates) are tagged too.
    auto split = WritePairs(ref, { { Aligned(30, 20, true), AlignmentResult() },
                                   { AlignmentResult(), Aligned(0, 20, false, 2) } },
                            Record("s/1", ref.gene.substr(30, 20)), Record("s/2", KmerUtils::ReverseComplement(ref.gene2.substr(0, 20))));
    ASSERT_EQ(split.size(), 2u);
    EXPECT_EQ(TagOf(split[0], "ZA"), "*");
    EXPECT_EQ(TagOf(split[1], "ZA"), "*");

    // Read back, a record keeps its ZA.
    Reader r(SamLine("a", 0).substr(0, SamLine("a", 0).size() - 1) + "\tZA:Z:2:1,7:0\n" + SamLine("b", 0));
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.sam1.m_alternatives, "2:1,7:0");
    EXPECT_NE(r.sam1.ToString().find("\tZA:Z:2:1,7:0"), std::string::npos);
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.sam1.m_alternatives, "");
    EXPECT_EQ(r.reader.PrimaryRecords(), 2u);
    EXPECT_EQ(r.reader.PrimaryRecordsWithoutAlternatives(), 1u);
}
