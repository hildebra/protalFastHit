// Unit tests for long reads: the chunks a read longer than seeding's 16-bit positions allow is
// seeded in, the segments of a read's hits, and a long read's SAM records as the profiler reads
// them back.
#include <gtest/gtest.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include "LongReads.h"
#include "Profiling/Profiler.h"
#include "ReadType.h"
#include "SequenceUtils/SeqReader.h"

using namespace protal;

TEST(ChunkRead, ReadsUpToTheLimitAreOneChunk) {
    auto chunks = ChunkRead(15000, kMaxLongReadChunk, 5000);
    ASSERT_EQ(chunks.size(), 1u);
    EXPECT_EQ(chunks[0].offset, 0u);
    EXPECT_EQ(chunks[0].length, 15000u);
    EXPECT_EQ(chunks[0].core.start, -kBeyondRead);  // genes the read starts or ends in are its own
    EXPECT_EQ(chunks[0].core.end, kBeyondRead);
    EXPECT_EQ(ChunkRead(kMaxLongReadChunk, kMaxLongReadChunk, 5000).size(), 1u);
}

TEST(ChunkRead, EveryGeneLiesInTheChunkThatOwnsIt) {
    size_t const max_chunk = 65000, overlap = 8000;
    for (size_t length : { 65001, 100000, 140000, 400003 }) {
        auto const chunks = ChunkRead(length, max_chunk, overlap);
        ASSERT_GT(chunks.size(), 1u) << length;
        EXPECT_EQ(chunks.front().core.start, -kBeyondRead);
        EXPECT_EQ(chunks.back().core.end, kBeyondRead);
        EXPECT_EQ(chunks.back().offset + chunks.back().length, length);
        for (size_t i = 0; i < chunks.size(); i++) {
            EXPECT_LE(chunks[i].length, max_chunk);
            if (i == 0) continue;
            EXPECT_EQ(chunks[i].core.start, chunks[i - 1].core.end) << "the cores tile the read";
            EXPECT_GE(chunks[i - 1].offset + chunks[i - 1].length - chunks[i].offset, overlap);
        }
        // Any stretch of up to `overlap` bases has one chunk whose core holds its centre, and lies in it
        // as far as it lies on the read (a gene may start before the read or end after it).
        for (int64_t start = -static_cast<int64_t>(overlap); start < static_cast<int64_t>(length); start += 997) {
            for (int64_t size : { int64_t{1}, int64_t{1500}, static_cast<int64_t>(overlap) }) {
                ReadInterval const gene{ start, start + size };
                if (gene.start >= static_cast<int64_t>(length) || gene.end <= 0) continue;
                auto owns = [&](ReadChunk const& c) { return gene.TwiceCentre() >= 2 * c.core.start && gene.TwiceCentre() < 2 * c.core.end; };
                ASSERT_EQ(std::count_if(chunks.begin(), chunks.end(), owns), 1) << length << " " << start;
                auto const& owner = *std::find_if(chunks.begin(), chunks.end(), owns);
                auto const on_read = gene.Within(static_cast<int64_t>(length));
                EXPECT_GE(on_read.start, static_cast<int64_t>(owner.offset));
                EXPECT_LE(on_read.end, static_cast<int64_t>(owner.offset + owner.length));
            }
        }
    }
}

namespace {
    LongReadCandidate Candidate(int64_t start, int64_t end) {
        return LongReadCandidate{ Anchor(1, 1, true), ReadChunk{}, ReadInterval{ start, end }, ReadInterval{} };
    }
}

TEST(LongReadSegments, HomologsShareASegmentAndNeighboursDoNot) {
    std::vector<LongReadCandidate> const candidates = {
            Candidate(1000, 2500),  // a gene
            Candidate(1050, 2450),  // its homolog in another taxon
            Candidate(2600, 3500),  // the next gene on the read
            Candidate(-300, 700),   // a gene the read starts in
            Candidate(2400, 3600),  // mostly on the next gene
    };
    auto const segments = LongReadSegmentsOf(candidates, 10000);
    ASSERT_EQ(segments.size(), 3u);
    EXPECT_EQ(segments[0], (std::vector<size_t>{ 0, 1 }));
    EXPECT_EQ(segments[1], (std::vector<size_t>{ 2, 4 }));
    EXPECT_EQ(segments[2], (std::vector<size_t>{ 3 }));
}

namespace {
    // A hit of `length` aligned bases on gene `geneid` of taxon `taxid`, `mismatches` of them counted as
    // mismatches (for the score only: the CIGAR is all M).
    LongReadHit Hit(uint32_t taxid, uint32_t geneid, size_t start, std::string seq, bool forward, size_t mismatches = 0) {
        LongReadHit hit;
        hit.alignment = AlignmentResult(0, taxid, geneid, static_cast<int32_t>(start), forward);
        auto& info = hit.alignment.GetAlignmentInfo();
        info.cigar = std::string(seq.size(), 'M');
        info.compressed_cigar = std::to_string(seq.size()) + "M";
        info.gene_alignment_start = static_cast<int>(start);
        info.alignment_length = seq.size();
        info.matches = seq.size() - mismatches;
        info.mismatches = mismatches;
        hit.qual = std::string(seq.size(), 'I');
        hit.seq = std::move(seq);
        return hit;
    }
}

TEST(LongReadSegments, RankedHitsCountOnceAndSetTheMapq) {
    LongReadSegment segment;
    segment.hits = { Hit(2, 1, 0, std::string(1000, 'A'), true, 10), Hit(1, 1, 0, std::string(1000, 'A'), true),
                     Hit(1, 1, 0, std::string(1000, 'A'), true) };
    RankLongReadSegment(segment);
    ASSERT_EQ(segment.hits.size(), 2u);
    EXPECT_EQ(segment.hits[0].alignment.Taxid(), 1u);
    EXPECT_EQ(segment.mapq, MAPQv2(2000, 990 * 2 - 10 * 3));
}

namespace {
    // A two-gene reference (taxid 1, genes 1 and 2) in a temporary directory, loaded.
    struct TinyReference {
        std::string gene = "ACGTTGCAAGGCTTACCGATGACTGAAACCGGTTTACGATCGGTAGCATG";   // 50 bp, gene 1
        std::string gene2 = "TTGACCAGTCAGGATCCATTGCAGGTACTTGACCGTAAGCTGCATTGACA";  // 50 bp, gene 2
        std::filesystem::path dir;
        std::unique_ptr<GenomeLoader> loader;

        TinyReference() {
            dir = std::filesystem::temp_directory_path() / ("protal_longread_test_" + std::to_string(::getpid()));
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

    // The records written; `genes`, if given, receives the genes the handler reported for them.
    std::string Write(TinyReference& ref, LongReadSegments segments, FastxRecord record, size_t max_out = 5,
                      std::vector<uint64_t>* genes = nullptr) {
        std::ostringstream os;
        SamStreamSink sink(os);
        {
            ProtalLongReadOutputHandler handler(sink, max_out, 1 << 16, *ref.loader);
            handler(segments, record);
        }  // the destructor flushes the buffer
        if (genes) *genes = sink.Genes();
        return os.str();
    }

    std::vector<std::vector<std::string>> Records(std::string const& sam) {
        std::vector<std::vector<std::string>> records;
        std::istringstream in(sam);
        std::string line;
        std::vector<std::string> tokens;
        while (std::getline(in, line)) {
            LineSplitter::Split(line, "\t", tokens);
            records.push_back(tokens);
        }
        return records;
    }

    FastxRecord Read(std::string id, std::string seq) {
        FastxRecord r;
        r.id = std::move(id);
        r.header = "@" + r.id;
        r.quality = std::string(seq.size(), 'I');
        r.sequence = std::move(seq);
        return r;
    }
}

TEST(LongReadOutputHandler, OnePrimaryAndSupplementaryRecordsHardClipped) {
    TinyReference ref;
    // Gene 1 forward at read positions 30-80, gene 2 reverse at 100-150, of a 175 bp read.
    std::string const read = std::string(30, 'C') + ref.gene + std::string(20, 'G') + KmerUtils::ReverseComplement(ref.gene2) + std::string(25, 'T');
    size_t const length = read.size();
    LongReadSegment first, second;
    first.hits = { Hit(1, 1, 0, ref.gene, true) };
    first.hits[0].clip_left = 30;
    first.hits[0].clip_right = length - 80;
    first.hits[0].aligned = { 30, 80 };
    first.mapq = 60;
    second.hits = { Hit(1, 2, 0, ref.gene2, false) };
    second.hits[0].clip_left = length - 150;  // on the reverse complement of the read
    second.hits[0].clip_right = 100;
    second.hits[0].aligned = { 100, 150 };
    second.mapq = 50;

    auto const sam = Write(ref, { second, first }, Read("movie/42/ccs", read));
    auto const records = Records(sam);
    ASSERT_EQ(records.size(), 2u);
    EXPECT_EQ(records[0][0], "movie/42/ccs");
    EXPECT_EQ(records[0][2], "1_1");
    EXPECT_EQ(std::stoul(records[0][1]), 0u);  // the primary alignment: the best segment, the first in read order of equals
    EXPECT_EQ(records[0][4], "60");
    EXPECT_EQ(records[0][5], "30H50M95H");
    EXPECT_EQ(records[0][9], ref.gene);
    EXPECT_EQ(records[1][2], "1_2");
    FLAG_t flag = std::stoul(records[1][1]);
    EXPECT_TRUE(Flag::IsSupplementaryAlignment(flag));
    EXPECT_TRUE(Flag::IsReverseComplement(flag));
    EXPECT_FALSE(Flag::IsNotPrimaryAlignment(flag));
    EXPECT_EQ(records[1][4], "50");
    EXPECT_EQ(records[1][5], "25H50M100H");
    EXPECT_EQ(records[1][9], ref.gene2);  // reference orientation

    // Read back, each segment is a read of its own.
    auto const path = (ref.dir / "long.sam").string();
    std::ofstream(path) << "@HD\tVN:1.6\n" << kSamReadTypeComment << "pb\n" << sam;
    profiler::Profiler profiler(*ref.loader);
    EXPECT_EQ(profiler.FromSam(path), "");
    ASSERT_EQ(profiler.m_pairs_unique.size(), 2u);
    EXPECT_EQ(profiler.m_pairs_unique[0].First().m_cigar, "50M");
    EXPECT_EQ(profiler.m_pairs_unique[1].First().m_rname, "1_2");
    std::ifstream in(path);
    auto const reads = ReadsOfSam(in);
    EXPECT_EQ(reads.declared, ReadType::PacBio);
}

TEST(LongReadOutputHandler, SecondaryHitsFollowTheirSegment) {
    TinyReference ref;
    // A read of Ns matches any position, so both hits are consistent.
    std::string const read(20, 'N');
    LongReadSegment segment;
    segment.hits = { Hit(1, 1, 5, read, true), Hit(1, 2, 3, read, true) };
    segment.mapq = 7;
    auto const records = Records(Write(ref, { segment }, Read("r", read)));
    ASSERT_EQ(records.size(), 2u);
    EXPECT_EQ(records[0][2], "1_1");
    EXPECT_EQ(records[0][4], "7");
    FLAG_t secondary = std::stoul(records[1][1]);
    EXPECT_TRUE(Flag::IsNotPrimaryAlignment(secondary));
    EXPECT_FALSE(Flag::IsSupplementaryAlignment(secondary));
    EXPECT_EQ(records[1][4], "0");
    EXPECT_EQ(Records(Write(ref, { segment }, Read("r", read), 1)).size(), 1u);

    // The genes of the records written, for the SAM header: both, or only the one written with -m 1.
    std::vector<uint64_t> genes;
    Write(ref, { segment }, Read("r", read), 5, &genes);
    EXPECT_EQ(genes, (std::vector<uint64_t>{ SamGeneKey(1, 1), SamGeneKey(1, 2) }));
    Write(ref, { segment }, Read("r", read), 1, &genes);
    EXPECT_EQ(genes, (std::vector<uint64_t>{ SamGeneKey(1, 1) }));
}

TEST(LongReadOutputHandler, SkipsOnlyTheInconsistentHit) {
    TinyReference ref;
    std::string const read = ref.gene.substr(0, 20);
    LongReadSegment segment;
    segment.hits = { Hit(1, 1, 3, read, true), Hit(1, 1, 0, read, true) };  // the first at the wrong position
    auto const records = Records(Write(ref, { segment }, Read("r", read)));
    ASSERT_EQ(records.size(), 1u);
    EXPECT_EQ(records[0][3], "1");
    EXPECT_EQ(std::stoul(records[0][1]), 0u);
}

namespace {
    // A segment whose best hit is on taxon `taxid` with MAPQ `mapq`.
    LongReadSegment Voting(uint32_t taxid, int mapq) {
        LongReadSegment segment;
        segment.hits = { Hit(taxid, 1, 0, std::string(100, 'A'), true) };
        segment.mapq = mapq;
        return segment;
    }
}

TEST(ReadConsensus, TaxaAreComparedOnTheGenesBothHit) {
    // Gene 1 fits taxon 1 better than taxon 2; gene 2 only taxon 2 (the database lacks taxon 1's gene 2, so the
    // read's gene aligns to the relative alone). Summed over the read, taxon 2 scores more (3500 against 2000), but
    // on the gene both hit it scores less: the read is taxon 1's.
    auto const read = ConsensusOfRead({ { { 1, 2000 }, { 2, 1750 } }, { { 2, 1750 } } });
    ASSERT_TRUE(read.has_value());
    EXPECT_EQ(read->taxid, 1u);
    EXPECT_EQ(read->parts, 1u);
    EXPECT_EQ(read->mapq, MAPQv2(2000, 1750));  // on taxon 1's gene, against taxon 2 there

    // A taxon without competitors on its genes: no MAPQ of the read (-1), its genes keep their own.
    auto const alone = ConsensusOfRead({ { { 1, 100 } }, { { 1, 120 } } });
    ASSERT_TRUE(alone.has_value());
    EXPECT_EQ(alone->mapq, -1);
    EXPECT_EQ(alone->parts, 2u);

    // Genes of three taxa that share no gene: none wins on enough genes.
    EXPECT_FALSE(ConsensusOfRead({ { { 1, 100 } }, { { 2, 100 } }, { { 3, 100 } } }, 2).has_value());
    EXPECT_FALSE(ConsensusOfRead({}).has_value());
    EXPECT_FALSE(ConsensusOfRead({ {}, {} }).has_value());

    // The read's genes: a taxon needs hits on half of them.
    auto segments = LongReadSegments{ Voting(1, 40), Voting(1, 30), Voting(2, 10) };
    auto const of_segments = ConsensusOfSegments(segments);
    ASSERT_TRUE(of_segments.has_value());
    EXPECT_EQ(of_segments->taxid, 1u);
}

TEST(ReadConsensus, EveryGeneTakesTheReadsTaxon) {
    ReadConsensus const read{ 1, 2, 30 };
    // The gene fits taxon 2 about as well as taxon 1 (1000 vs 998 matches): settled for taxon 1.
    LongReadSegment ambiguous;
    ambiguous.hits = { Hit(2, 7, 0, std::string(1000, 'A'), true), Hit(1, 7, 0, std::string(1000, 'A'), true, 2) };
    RankLongReadSegment(ambiguous);
    ASSERT_LT(ambiguous.mapq, kConfidentMapq);
    EXPECT_TRUE(SettleByRead(ambiguous, read));
    EXPECT_EQ(ambiguous.hits.front().alignment.Taxid(), 1u);
    EXPECT_EQ(ambiguous.mapq, 30);
    EXPECT_TRUE(ambiguous.by_read);

    // A gene that fits taxon 2 clearly better keeps it (a chimeric read, or a homolog of a gene of taxon 1 that lies
    // elsewhere on the read), with MAPQ 0: no evidence for taxon 2.
    LongReadSegment other;
    other.hits = { Hit(2, 7, 0, std::string(1000, 'A'), true), Hit(1, 7, 0, std::string(1000, 'A'), true, 100) };
    RankLongReadSegment(other);
    ASSERT_GE(other.mapq, kConfidentMapq);
    EXPECT_TRUE(SettleByRead(other, read));
    EXPECT_EQ(other.hits.front().alignment.Taxid(), 2u);
    EXPECT_EQ(other.mapq, 0);
    EXPECT_TRUE(other.inconsistent);
    EXPECT_FALSE(other.by_read);

    // A gene of the read's taxon keeps a higher MAPQ of its own, and a lower one is raised.
    auto confident = Voting(1, 40);
    EXPECT_FALSE(SettleByRead(confident, read));
    EXPECT_EQ(confident.mapq, 40);
    EXPECT_FALSE(confident.by_read);
    auto weak = Voting(1, 2);
    EXPECT_TRUE(SettleByRead(weak, read));
    EXPECT_EQ(weak.mapq, 30);
    EXPECT_TRUE(weak.by_read);
    // Without competitors on the read (MAPQ -1) a gene keeps its own.
    auto own = Voting(1, 5);
    EXPECT_FALSE(SettleByRead(own, ReadConsensus{ 1, 2, -1 }));
    EXPECT_EQ(own.mapq, 5);

    // A gene without a hit of the read's taxon is no evidence for its own: MAPQ 0, inconsistent.
    auto alone = Voting(2, 60);
    EXPECT_TRUE(SettleByRead(alone, read));
    EXPECT_EQ(alone.mapq, 0);
    EXPECT_TRUE(alone.inconsistent);
    EXPECT_FALSE(alone.by_read);
    EXPECT_EQ(alone.hits.front().alignment.Taxid(), 2u);
}

TEST(LongReadOutputHandler, TagsGenesSettledByTheirRead) {
    TinyReference ref;
    std::string const read(20, 'N');
    LongReadSegment settled;
    settled.hits = { Hit(1, 1, 5, read, true), Hit(1, 2, 3, read, true) };
    settled.mapq = 30;
    settled.by_read = true;
    auto const records = Records(Write(ref, { settled }, Read("r", read)));
    ASSERT_EQ(records.size(), 2u);
    EXPECT_EQ(records[0].back(), "ZR:i:1");
    EXPECT_NE(records[1].back(), "ZR:i:1") << "only the best hit";

    // A gene without a hit of its read's taxon: ZR:i:2, MAPQ 0.
    LongReadSegment inconsistent;
    inconsistent.hits = { Hit(1, 1, 5, read, true) };
    inconsistent.mapq = 0;
    inconsistent.inconsistent = true;
    auto const tagged = Records(Write(ref, { inconsistent }, Read("r", read)));
    ASSERT_EQ(tagged.size(), 1u);
    EXPECT_EQ(tagged[0].back(), "ZR:i:2");
    EXPECT_EQ(tagged[0][4], "0");
}

TEST(SeqReader, ReadsWithoutQualitiesGetTheReadTypesQuality) {
    std::istringstream fasta(">r1\nACGT\n>r2\nAC\nGT\n");
    SeqReaderSE reader(fasta, FastaQualityChar(ReadType::PacBio));
    FastxRecord record;
    ASSERT_TRUE(reader(record));
    EXPECT_EQ(record.quality, "????");  // Q30
    ASSERT_TRUE(reader(record));
    EXPECT_EQ(record.sequence, "ACGT");
    EXPECT_EQ(record.quality, "????");
    EXPECT_FALSE(reader(record));

    std::istringstream fastq("@r1\nACGT\n+\nIII#\n");
    SeqReaderSE keeps(fastq, FastaQualityChar(ReadType::PacBio));
    ASSERT_TRUE(keeps(record));
    EXPECT_EQ(record.quality, "III#") << "FASTQ keeps its qualities";

    std::istringstream mate1(">p/1\nAAAA\n"), mate2(">p/2\nCCC\n");
    SeqReaderPE pairs(mate1, mate2, FastaQualityChar(ReadType::Paired));
    FastxRecord record2;
    ASSERT_TRUE(pairs(record, record2));
    EXPECT_EQ(record.quality, "????");
    EXPECT_EQ(record2.quality, "???");
}
