// Unit tests for parsing tab-separated input: the line splitter, sample map rows, and SAM records
// as the profiler reads them (protal's own and other aligners').
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include "Options.h"
#include "Profiling/Profiler.h"
#include "IO/SamHandler.h"
#include "Utilities/LineSplitter.h"

using namespace protal;
namespace fs = std::filesystem;

namespace {
    std::vector<std::string> Split(std::string const& line, std::string const& delimiter = "\t") {
        std::vector<std::string> tokens;
        LineSplitter::Split(line, delimiter, tokens);
        return tokens;
    }

    using Tokens = std::vector<std::string>;
}

TEST(LineSplitter, KeepsEmptyFields) {
    EXPECT_EQ(Split("a\tb\tc"), (Tokens{ "a", "b", "c" }));
    EXPECT_EQ(Split("a\t\tb"), (Tokens{ "a", "", "b" }));
    EXPECT_EQ(Split("a\t"), (Tokens{ "a", "" }));
    EXPECT_EQ(Split("\ta"), (Tokens{ "", "a" }));
    EXPECT_EQ(Split("\t"), (Tokens{ "", "" }));
    EXPECT_EQ(Split("a"), (Tokens{ "a" }));
    EXPECT_TRUE(Split("").empty());
    EXPECT_EQ(Split("a::b::::c", "::"), (Tokens{ "a", "b", "", "c" }));
    EXPECT_EQ(Split("s1.fq,,s3.fq", ","), (Tokens{ "s1.fq", "", "s3.fq" }));
}

namespace {
    struct ScratchDir {
        fs::path path;
        ScratchDir() {
            path = fs::temp_directory_path() / ("protal parsing test " + std::to_string(::getpid()));
            fs::create_directories(path);
        }
        ~ScratchDir() { fs::remove_all(path); }

        std::string Write(std::string const& name, std::string const& content) const {
            auto file = path / name;
            std::ofstream(file, std::ios::binary) << content;
            return file.string();
        }
    };

    struct MapLists {
        std::string output_dir, strain_dir, misc_dir;
        std::vector<std::string> prefixes, firsts, seconds, sams, profiles, names, truths, read_types;

        bool Load(std::string const& path) {
            return Options::LoadFromMap(path, output_dir, strain_dir, misc_dir, prefixes, firsts, seconds, sams, profiles, names, truths, read_types);
        }
    };

    std::string MapHeader(fs::path const& out) {
        return "#OUTPUT_DIR\t" + out.string() + "\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tPROFILE\n";
    }
}

TEST(SampleMap, ReadsRowsWithCRLF) {
    ScratchDir dir;
    std::string map = MapHeader(dir.path / "out") + "s1\ts1\ta_1.fq\ta_2.fq\ts1.profile\r\n\ns2\ts2\tb_1.fq\tb_2.fq\ts2.profile\n";
    MapLists lists;
    ASSERT_TRUE(lists.Load(dir.Write("samples.map", map)));
    EXPECT_EQ(lists.names, (Tokens{ "s1", "s2" }));
    ASSERT_EQ(lists.profiles.size(), 2u);
    EXPECT_TRUE(lists.profiles[1].ends_with("s2.profile"));
    EXPECT_TRUE(lists.seconds[0].ends_with("a_2.fq"));
}

TEST(SampleMap, RejectsRowsWithMissingOrEmptyCells) {
    ScratchDir dir;
    auto header = MapHeader(dir.path / "out");
    testing::internal::CaptureStderr();
    MapLists short_row;
    EXPECT_FALSE(short_row.Load(dir.Write("short.map", header + "s1\ts1\ta_1.fq\ta_2.fq\ts1.profile\ns2\ts2\tb_1.fq\tb_2.fq\n")));
    MapLists empty_cell;
    EXPECT_FALSE(empty_cell.Load(dir.Write("empty.map", header + "s1\ts1\t\ta_2.fq\ts1.profile\n")));
    auto log = testing::internal::GetCapturedStderr();
    EXPECT_NE(log.find("Line 4: no value in column 5 (PROFILE)"), std::string::npos) << log;
    EXPECT_NE(log.find("Line 3: no value in column 3 (FIRST)"), std::string::npos) << log;
}

TEST(SampleMap, SingleEndSamplesHaveNoSecondFile) {
    ScratchDir dir;
    auto const out = dir.path / "out";
    // '-' in SECOND marks a single-end sample among paired-end ones.
    MapLists mixed;
    ASSERT_TRUE(mixed.Load(dir.Write("mixed.map", "#OUTPUT_DIR\t" + out.string() + "\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\n"
                                                  "p\tp\tp_1.fq\tp_2.fq\ns\ts\ts.fq\t-\n")));
    ASSERT_EQ(mixed.seconds.size(), 2u);
    EXPECT_TRUE(mixed.seconds[0].ends_with("p_2.fq"));
    EXPECT_EQ(mixed.seconds[1], "");
    EXPECT_TRUE(mixed.firsts[1].ends_with("s.fq"));

    // Without a SECOND column, all samples are single-end.
    MapLists single;
    ASSERT_TRUE(single.Load(dir.Write("single.map", "#OUTPUT_DIR\t" + out.string() + "\n#SAMPLEID\tPREFIX\tFIRST\n"
                                                    "a\ta\ta.fq.gz\nb\tb\tb.fq.gz\n")));
    EXPECT_EQ(single.seconds, (Tokens{ "", "" }));
    EXPECT_EQ(single.firsts.size(), 2u);
    EXPECT_EQ(single.sams.size(), 2u);
    EXPECT_TRUE(single.read_types.empty());  // --read_type applies

    // READ_TYPE names each sample's reads.
    MapLists typed;
    ASSERT_TRUE(typed.Load(dir.Write("typed.map", "#OUTPUT_DIR\t" + out.string() + "\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tREAD_TYPE\n"
                                                  "p\tp\tp_1.fq\tp_2.fq\tpe\nl\tl\tl.fq.gz\t-\tpb\n")));
    EXPECT_EQ(typed.read_types, (Tokens{ "pe", "pb" }));
    EXPECT_EQ(typed.seconds[1], "");
}

TEST(ReadFileStem, DropsReadAndCompressionExtensions) {
    EXPECT_EQ(Utils::ReadFileStem("sample1.fq.gz"), "sample1");
    EXPECT_EQ(Utils::ReadFileStem("sample1_R1.fastq"), "sample1_R1");
    EXPECT_EQ(Utils::ReadFileStem("reads.fa.zst"), "reads");
    EXPECT_EQ(Utils::ReadFileStem("reads.gz"), "reads");
    EXPECT_EQ(Utils::ReadFileStem("reads.txt"), "reads.txt");
    EXPECT_EQ(Utils::ReadFileStem(".fq"), ".fq");
}

TEST(Options, EachReadTypeHasItsModel) {
    ScratchDir dir;
    auto model = [&](std::string const& paired, std::string const& own, ReadType type) {
        OptionsData d;
        d.database_path = dir.path.string();
        d.model = paired;
        (type == ReadType::PacBio ? d.model_pb : type == ReadType::ONT ? d.model_ont : d.model_se) = own;
        return fs::path(Options(d).ModelDbFile(type).Path()).filename().string();
    };
    EXPECT_EQ(model("", "", ReadType::Paired), "model_pe.xml");
    dir.Write("model.xml", "");  // a database from before read types
    EXPECT_EQ(model("", "", ReadType::Paired), "model.xml");
    dir.Write("model_pe.xml", "");
    EXPECT_EQ(model("", "", ReadType::Paired), "model_pe.xml");
    EXPECT_EQ(model("", "", ReadType::Single), "model_se.xml");
    EXPECT_EQ(model("", "", ReadType::PacBio), "model_PB.xml");
    // --model replaces all, unless --model_se or --model_pb is given for their reads.
    EXPECT_EQ(model("other", "", ReadType::Paired), "other.xml");
    EXPECT_EQ(model("other", "", ReadType::Single), "other.xml");
    EXPECT_EQ(model("other", "", ReadType::PacBio), "other.xml");
    EXPECT_EQ(model("other", "se2.xml", ReadType::Single), "se2.xml");
    EXPECT_EQ(model("other", "se2.xml", ReadType::Paired), "other.xml");
    EXPECT_EQ(model("other", "hifi", ReadType::PacBio), "hifi.xml");
    auto const existing = dir.Write("elsewhere.xml", "");
    EXPECT_EQ(model("", existing, ReadType::Single), "elsewhere.xml");
    EXPECT_EQ(model("", "", ReadType::ONT), "model_ONT.xml");
    EXPECT_EQ(model("other", "", ReadType::ONT), "other.xml");
    EXPECT_EQ(model("other", "r10", ReadType::ONT), "r10.xml");
}

TEST(Options, OntReadsHaveTheirOwnDefaults) {
    OptionsData d;
    Options defaults(d);
    EXPECT_DOUBLE_EQ(defaults.GetMaxScoreAni(ReadType::ONT), 0.85);
    EXPECT_DOUBLE_EQ(defaults.GetSNPMinAF(ReadType::ONT), 0.2);
    EXPECT_DOUBLE_EQ(defaults.GetMaxScoreAni(ReadType::PacBio), DEFAULT_MAX_SCORE_ANI);
    EXPECT_DOUBLE_EQ(defaults.GetSNPMinAF(ReadType::Paired), DEFAULT_MIN_SNP_AF);
    EXPECT_EQ(FastaQualityChar(ReadType::ONT), '3');  // Q18
    EXPECT_EQ(FastaQualityChar(ReadType::Single), '?');  // Q30
    // Given, the options apply to all read types.
    d.max_score_ani = 0.95;
    d.max_score_ani_given = true;
    d.snp_min_af = 0.1;
    d.snp_min_af_given = true;
    Options given(d);
    EXPECT_DOUBLE_EQ(given.GetMaxScoreAni(ReadType::ONT), 0.95);
    EXPECT_DOUBLE_EQ(given.GetSNPMinAF(ReadType::ONT), 0.1);
}

namespace {
    std::string Record(std::string const& qname, int flag, std::string const& rname, std::string const& cigar,
                       std::string const& seq, std::string const& tags = "\tZU:i:1\tZT:i:0", int pos = 1, int mapq = 60) {
        return qname + '\t' + std::to_string(flag) + '\t' + rname + '\t' + std::to_string(pos) + '\t' + std::to_string(mapq) + '\t' +
               cigar + "\t*\t0\t0\t" + seq + '\t' + (seq == "*" ? std::string("*") : std::string(seq.size(), 'I')) + tags + '\n';
    }

    struct Reader {
        std::istringstream in;
        SamReader reader{ in };
        SamEntry sam1, sam2;
        bool has1 = false, has2 = false;

        explicit Reader(std::string text) : in(std::move(text)) {}
        bool Next() { return reader.Next(sam1, sam2, has1, has2); }
        SamEntry const& Any() const { return has1 ? sam1 : sam2; }
    };
}

TEST(SamReader, SkipsHeadersBlankLinesAndUnusableRecords) {
    Reader r("@HD\tVN:1.6\n@SQ\tSN:1_1\tLN:50\n\n" +
             Record("a", 0, "1_1", "4M", "ACGT") + "\r\n" +
             Record("unmapped", 4, "1_1", "4M", "ACGT") +
             Record("star", 0, "*", "*", "ACGT") +
             Record("chr", 0, "chr1", "4M", "ACGT") +
             Record("noseq", 256, "1_1", "4M", "*") +
             Record("spliced", 0, "1_1", "2M5N2M", "ACGT") +
             Record("short", 0, "1_1", "3M", "ACGT") +
             "\n" + Record("b", 0, "2_7", "4M", "ACGT"));
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.Any().m_qname, "a");
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.Any().m_qname, "b");
    EXPECT_FALSE(r.Next());

    auto const& skipped = r.reader.Skipped();
    EXPECT_EQ(r.reader.Records(), 2u);
    EXPECT_EQ(skipped.at("unmapped"), 2u);
    EXPECT_EQ(skipped.at("reference is not a protal gene (<taxid>_<gene id>)"), 1u);
    EXPECT_EQ(skipped.at("no sequence"), 1u);
    EXPECT_EQ(skipped.at("CIGAR is missing, malformed, has N/P ops or does not match SEQ"), 2u);
}

TEST(SamReader, LastRecordWithoutNewlineIsRead) {
    auto text = Record("a", 0, "1_1", "4M", "ACGT");
    text.pop_back();
    Reader r("@HD\tVN:1.6\n" + text);
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.Any().m_qname, "a");
    EXPECT_FALSE(r.Next());
}

TEST(SamReader, NormalizesSequenceMatchAndHardClips) {
    Reader r(Record("a", 0, "1_1", "5H2=1X1=3H", "ACGT") + Record("b", 0, "1_1", "1S2M1M", "ACGT"));
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.Any().m_cigar, "2M1X1M");
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.Any().m_cigar, "1S3M");

    std::string cigar = "4M";
    EXPECT_FALSE(NormalizeCigar(cigar, 5));
    cigar = "0M4M";
    EXPECT_FALSE(NormalizeCigar(cigar, 4));
    cigar = "4";
    EXPECT_FALSE(NormalizeCigar(cigar, 4));
    cigar = "2M2P2M";
    EXPECT_FALSE(NormalizeCigar(cigar, 4));
}

TEST(SamReader, ReadsTagsByName) {
    Reader r(Record("protal", 0, "1_1", "4M", "ACGT", "\tZU:i:7\tZT:i:2") +
             Record("reordered", 0, "1_1", "4M", "ACGT", "\tNM:i:3\tRG:Z:x\tZT:i:4\tZU:i:9") +
             Record("bowtie2", 0, "1_1", "4M", "ACGT", "\tAS:i:-5\tXS:i:-12") +
             Record("none", 0, "1_1", "4M", "ACGT", ""));
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.Any().m_uniques, 7);
    EXPECT_EQ(r.Any().m_uniques_two, 2);
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.Any().m_uniques, 9);
    EXPECT_EQ(r.Any().m_uniques_two, 4);
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.Any().m_uniques, 0);
    ASSERT_TRUE(r.Next());
    EXPECT_EQ(r.Any().m_uniques, 0);
    EXPECT_FALSE(r.Next());
    EXPECT_EQ(r.reader.RecordsWithoutTags(), 2u);
}

TEST(SamReader, MalformedLinesThrowWithTheLineNumber) {
    auto expect_error = [](std::string const& text, std::string const& message) {
        Reader r(text);
        try {
            while (r.Next()) {}
            ADD_FAILURE() << "no error for: " << text;
        } catch (SamFormatError const& e) {
            EXPECT_NE(std::string(e.what()).find(message), std::string::npos) << e.what();
        }
    };
    expect_error("@HD\tVN:1.6\n" + Record("a", 0, "1_1", "4M", "ACGT") + "b\t0\t1_1\t1\n", "line 3: expected at least 11 tab-separated fields, found 4");
    expect_error("a\tx\t1_1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n", "line 1: FLAG is not a 16-bit integer: x");
    expect_error("a\t0\t1_1\t1\t300\t4M\t*\t0\t0\tACGT\tIIII\n", "MAPQ is not between 0 and 255");
    expect_error("a\t0\t1_1\t-1\t60\t4M\t*\t0\t0\tACGT\tIIII\n", "POS is not a position");
}

namespace {
    // A one-gene reference (taxid 1, gene 1).
    struct TinyReference {
        std::string gene = "ACGTTGCAAGGCTTACCGATGACTGAAACCGGTTTACGATCGGTAGCATG";  // 50 bp
        ScratchDir dir;
        std::unique_ptr<GenomeLoader> loader;

        TinyReference() {
            std::string header = ">1_1\n";
            auto fna = dir.Write("reference.fna", header + gene + '\n');
            auto map = dir.Write("reference.map", "1\t1\t" + std::to_string(header.size()) + '\t' + std::to_string(header.size() + gene.size()) + '\n');
            loader = std::make_unique<GenomeLoader>(fna, map);
            loader->LoadAllGenomes();
        }
    };

    constexpr int kPaired = 0x1, kBothAlign = 0x2, kRead1 = 0x40, kRead2 = 0x80, kSecondary = 0x100;
}

TEST(FromSam, GroupsCandidatesAndKeepsTheLastGroup) {
    TinyReference ref;
    auto seq = ref.gene.substr(0, 20);
    int const r1 = kPaired | kBothAlign | kRead1, r2 = kPaired | kBothAlign | kRead2;
    auto sam = ref.dir.Write("sample.sam", "@HD\tVN:1.6\n" +
            Record("u", r1, "1_1", "20M", seq) + Record("u", r2, "1_1", "20M", seq) +
            // The primary alignment is the second candidate here.
            Record("m", r1 | kSecondary, "1_1", "20M", seq, "\tZU:i:1", 1, 0) + Record("m", r2 | kSecondary, "1_1", "20M", seq, "\tZU:i:1", 1, 0) +
            Record("m", r1, "1_1", "20M", seq, "\tZU:i:1", 5) + Record("m", r2, "1_1", "20M", seq, "\tZU:i:1", 5) +
            // The file ends in a multi-mapped read.
            Record("z", 0, "1_1", "20M", seq) + Record("z", kSecondary, "1_1", "20M", seq, "\tZU:i:1", 9, 0));

    profiler::Profiler profiler(*ref.loader);
    EXPECT_EQ(profiler.FromSam(sam), "");
    ASSERT_EQ(profiler.m_pairs_unique.size(), 1u);
    EXPECT_EQ(profiler.m_pairs_unique[0].First().m_qname, "u");
    EXPECT_TRUE(profiler.m_pairs_unique[0].HasSecond());

    ASSERT_EQ(profiler.m_pairs_nonunique.size(), 2u);
    ASSERT_EQ(profiler.m_pairs_nonunique_best.size(), 2u);
    EXPECT_EQ(profiler.m_pairs_nonunique[0].size(), 2u);
    EXPECT_EQ(profiler.m_pairs_nonunique_best[0].First().m_pos, 5u);
    EXPECT_EQ(profiler.m_pairs_nonunique[1].size(), 2u);
    EXPECT_EQ(profiler.m_pairs_nonunique_best[1].Any().m_qname, "z");
    EXPECT_EQ(profiler.m_pairs_nonunique_best[1].Any().m_pos, 1u);
}

TEST(FromSam, EmptyHeaderOnlyAndBrokenFiles) {
    TinyReference ref;
    profiler::Profiler profiler(*ref.loader);

    testing::internal::CaptureStderr();
    EXPECT_EQ(profiler.FromSam(ref.dir.Write("header.sam", "@HD\tVN:1.6\n@SQ\tSN:1_1\tLN:50\n")), "");
    EXPECT_FALSE(profiler.HasReads());
    auto log = testing::internal::GetCapturedStderr();
    EXPECT_NE(log.find("contains no usable alignments"), std::string::npos);

    EXPECT_NE(profiler.FromSam(ref.dir.Write("empty.sam", "")).find("the file is empty"), std::string::npos);
    auto broken = profiler.FromSam(ref.dir.Write("broken.sam", "@HD\tVN:1.6\n" + Record("a", 0, "1_1", "4M", "ACGT") + "b\t0\t1_1"));
    EXPECT_NE(broken.find("line 3"), std::string::npos) << broken;
}

TEST(FromSam, RejectsASamAlignedAgainstAnotherDatabase) {
    TinyReference ref;
    profiler::Profiler profiler(*ref.loader);
    auto record = Record("a", 0, "1_1", "20M", ref.gene.substr(0, 20));
    auto with_header = [&](std::string const& sq) {
        return profiler.FromSam(ref.dir.Write("other.sam", "@HD\tVN:1.6\n" + sq + record));
    };
    EXPECT_EQ(with_header("@SQ\tSN:1_1\tLN:50\n"), "");
    EXPECT_EQ(with_header(""), "");  // no header to check against
    EXPECT_EQ(with_header("@SQ\tSN:chr1\tLN:1000\n"), "");  // not a protal gene: its records are skipped

    auto length = with_header("@SQ\tSN:1_1\tLN:60\n");
    EXPECT_NE(length.find("gene 1_1 is 60 bp in the SAM header (@SQ) but 50 bp in the database"), std::string::npos) << length;
    auto missing = with_header("@SQ\tSN:1_1\tLN:50\n@SQ\tSN:7_3\tLN:50\n");
    EXPECT_NE(missing.find("gene 7_3 of the SAM header (@SQ) is not in the database"), std::string::npos) << missing;
}

TEST(FromSam, ATruncatedGzipFileIsAnError) {
    TinyReference ref;
    std::string sam = "@HD\tVN:1.6\n@SQ\tSN:1_1\tLN:50\n";
    for (int i = 0; i < 20000; i++) sam += Record("r" + std::to_string(i), 0, "1_1", "20M", ref.gene.substr(i % 30, 20), "\tZU:i:1", 1 + i % 30);
    auto path = (std::filesystem::path(ref.dir.path) / "sample.sam.gz").string();
    {
        ogzstream os(path.c_str());
        os << sam;
    }
    profiler::Profiler profiler(*ref.loader);
    EXPECT_EQ(profiler.FromSam(path), "");
    EXPECT_EQ(profiler.m_pairs_unique.size(), 20000u);

    std::filesystem::resize_file(path, std::filesystem::file_size(path) / 2);
    auto error = profiler.FromSam(path);
    EXPECT_NE(error.find("the file is truncated or corrupt"), std::string::npos) << error;
}

TEST(FromSam, ReadsTheCompressedSamsProtalWrites) {
    TinyReference ref;
    std::string const header = "@HD\tVN:1.6\n@SQ\tSN:1_1\tLN:50\n";
    std::string records;
    for (int i = 0; i < 20000; i++) records += Record("r" + std::to_string(i), 0, "1_1", "20M", ref.gene.substr(i % 30, 20), "\tZU:i:1", 1 + i % 30);
    for (auto const& name : { "sample.sam.gz", "sample.sam.zst" }) {
        SCOPED_TRACE(name);
        auto const path = (ref.dir.path / name).string();
        {
            SamOutput out(path, SamCompressionOf(path));
            std::vector<uint64_t> genes = { SamGeneKey(1, 1) };
            out.Write(records.data(), records.size(), genes);
            ASSERT_TRUE(out.Finish(header)) << out.Error();
        }
        profiler::Profiler profiler(*ref.loader);
        EXPECT_EQ(profiler.FromSam(path), "");
        EXPECT_EQ(profiler.m_pairs_unique.size(), 20000u);

        // Cut after its last block or frame, the file still decompresses, but without its end
        // marker (BGZF's end-of-file block, the seek table): an error, not a sample with fewer reads.
        uint64_t cut = fs::file_size(path) - 28;
        if (SamCompressionOf(path) == SamCompression::Zstd) {
            std::string error;
            auto const table = zstd::ReadSeekTable(path, error);
            ASSERT_TRUE(table.has_value()) << error;
            cut = table->frames.back().compressed_offset + table->frames.back().compressed_size;
        }
        fs::resize_file(path, cut);
        auto const error = profiler.FromSam(path);
        EXPECT_NE(error.find("the file is truncated or corrupt"), std::string::npos) << error;
    }
}

TEST(Options, ReadTypesOfSamples) {
    ScratchDir dir;
    // From the read files: no second file means single-end reads.
    OptionsData reads;
    // PacBio reads come as READ_TYPE (or --read_type) pb.
    reads.first_list = { "p_1.fq", "s.fq", "l.fq" };
    reads.second_list = { "p_2.fq", "", "" };
    reads.prefix_list = { "p", "s", "l" };
    reads.read_type_list = { "", "", "pb" };
    reads.range = { 0, 1 };
    Options from_reads(reads);
    std::vector<std::string> warnings;
    from_reads.ResolveReadTypes(warnings);
    EXPECT_EQ(from_reads.GetReadType(0), ReadType::Paired);
    EXPECT_EQ(from_reads.GetReadType(1), ReadType::Single);
    EXPECT_TRUE(from_reads.AnySample(ReadType::Single) && from_reads.AnySample(ReadType::Paired));
    EXPECT_EQ(from_reads.GetReadType(2), ReadType::PacBio);
    EXPECT_FALSE(from_reads.AnySample(ReadType::PacBio));  // not in the range

    // With --profile_only, from the SAM: the kind its header names, else unpaired records are
    // single-end reads; a SAM without alignments counts as paired-end.
    auto seq = std::string("ACGTACGTACGTACGTACGT");
    OptionsData sams;
    sams.profile_only = true;
    sams.sam_list = { dir.Write("pe.sam", "@HD\tVN:1.6\n" + Record("a", kPaired | kBothAlign | kRead1, "1_1", "20M", seq) +
                                          Record("a", kPaired | kBothAlign | kRead2, "1_1", "20M", seq)),
                      dir.Write("se.sam", "@HD\tVN:1.6\n" + Record("b", 16, "1_1", "20M", seq)),
                      dir.Write("empty.sam", "@HD\tVN:1.6\n"),
                      dir.Write("long.sam", "@HD\tVN:1.6\n" + kSamReadTypeComment + "pb\n" + Record("c", 0x800, "1_1", "5H20M", seq)) };
    sams.prefix_list = { "pe", "se", "empty", "long" };
    sams.range = { 1 };
    Options from_sams(sams);
    from_sams.ResolveReadTypes(warnings);
    EXPECT_EQ(from_sams.GetReadType(0), ReadType::Paired);
    EXPECT_EQ(from_sams.GetReadType(1), ReadType::Single);
    EXPECT_EQ(from_sams.GetReadType(2), ReadType::Paired);
    EXPECT_EQ(from_sams.GetReadType(3), ReadType::PacBio);
    EXPECT_TRUE(from_sams.AnySample(ReadType::Single));
    EXPECT_FALSE(from_sams.AnySample(ReadType::Paired));  // only the single-end sample is in the range
    EXPECT_TRUE(warnings.empty());

    // A read type given for a SAM wins over its own, with a warning where they differ.
    sams.read_type_list = { "se", "se", "", "pb" };
    Options given(sams);
    given.ResolveReadTypes(warnings);
    EXPECT_EQ(given.GetReadType(0), ReadType::Single);
    EXPECT_EQ(given.GetReadType(1), ReadType::Single);
    EXPECT_EQ(given.GetReadType(2), ReadType::Paired);
    EXPECT_EQ(given.GetReadType(3), ReadType::PacBio);
    ASSERT_EQ(warnings.size(), 1u);
    EXPECT_NE(warnings[0].find("pe.sam holds paired-end reads; profiled as single-end reads"), std::string::npos) << warnings[0];
}

TEST(MicrobialProfile, RejectsRecordsOutsideTheDatabase) {
    TinyReference ref;
    profiler::MicrobialProfile profile(*ref.loader);
    auto sam_at = [&](std::string rname, int pos, std::string const& cigar, std::string const& seq) {
        SamEntry sam;
        sam.m_qname = "r";
        sam.m_flag = 0;
        sam.m_rname = std::move(rname);
        sam.m_pos = pos;
        sam.m_mapq = 60;
        sam.m_cigar = cigar;
        sam.m_seq = seq;
        sam.m_qual = std::string(seq.size(), 'I');
        return sam;
    };
    auto foreign_taxon = sam_at("9_1", 1, "20M", ref.gene.substr(0, 20));
    auto foreign_gene = sam_at("1_2", 1, "20M", ref.gene.substr(0, 20));
    auto past_end = sam_at("1_1", 41, "20M", ref.gene.substr(30, 20));
    EXPECT_FALSE(profile.AddSam(9, 1, foreign_taxon, 1.0));
    EXPECT_FALSE(profile.AddSam(1, 2, foreign_gene, 1.0));
    EXPECT_FALSE(profile.AddSam(1, 1, past_end, 1.0));
    EXPECT_TRUE(profile.GetTaxa().empty());

    auto inside = sam_at("1_1", 31, "20M", ref.gene.substr(30, 20));
    EXPECT_TRUE(profile.AddSam(1, 1, inside, 1.0));
    ASSERT_EQ(profile.GetTaxa().size(), 1u);
    EXPECT_EQ(profile.GetTaxa().at(1).TotalHits(), 1u);
}

TEST(MicrobialProfile, CountsAlternativesLowMapqAndLinkedReads) {
    TinyReference ref;
    auto sam_at = [&](int pos, int mapq, std::string alternatives) {
        SamEntry sam;
        sam.m_qname = "r";
        sam.m_flag = 0;
        sam.m_rname = "1_1";
        sam.m_pos = pos;
        sam.m_mapq = mapq;
        sam.m_cigar = "10M";
        sam.m_seq = ref.gene.substr(pos - 1, 10);
        sam.m_qual = std::string(10, 'I');
        sam.m_alternatives = std::move(alternatives);
        return sam;
    };
    // Taxa 1 and 2 are of genus 10, taxon 3 of genus 11.
    auto genera = std::make_shared<std::vector<uint32_t>>(std::vector<uint32_t>{ 0, 10, 10, 11 });
    auto fill = [&](profiler::MicrobialProfile& profile) {
        auto add = [&](SamEntry const& sam, int read) {
            profile.NoteRecord(1, sam);
            EXPECT_TRUE(profile.AddSam(1, 1, sam, 1.0, true, read, true, static_cast<size_t>(read)));
        };
        // Read 0, both mates here: one fits a congener within an edit, the other a species of another genus.
        add(sam_at(1, 60, "2:1"), 0);
        add(sam_at(21, 60, "3:0,2:4"), 0);
        // Read 1: an alternative two edits worse does not fit as well; its MAPQ is low.
        add(sam_at(5, 3, "2:2"), 1);
        // Read 2: no alternative.
        add(sam_at(31, 60, "*"), 2);
        // Read 3 fits a congener exactly (MAPQ 0): the profiler's MAPQ filter leaves it out of the taxon's hits,
        // but it counts for the evidence.
        profile.NoteRecord(1, sam_at(11, 0, "2:0"));
        profile.ApplyRecordEvidence();
        return profile.GetTaxa().at(1);
    };
    profiler::MicrobialProfile profile(*ref.loader);
    profile.SetGenera(genera);
    auto const& taxon = fill(profile);
    EXPECT_EQ(taxon.TotalHits(), 4u);
    EXPECT_DOUBLE_EQ(taxon.CongenerFitShare(), 0.4);
    EXPECT_DOUBLE_EQ(taxon.OtherGenusFitShare(), 0.2);
    EXPECT_DOUBLE_EQ(taxon.LowMapqShare(), 0.4);
    EXPECT_DOUBLE_EQ(taxon.LinkedShare(), 1.0 / 3);
    // As the model and the training dump get them.
    std::map<std::string, double> features;
    for (auto const& [name, value] : profiler::TaxonFeatures(taxon)) features[name] = value;
    EXPECT_DOUBLE_EQ(features.at("congener_fit_share"), 0.4);
    EXPECT_DOUBLE_EQ(features.at("other_genus_fit_share"), 0.2);
    EXPECT_DOUBLE_EQ(features.at("low_mapq_share"), 0.4);
    EXPECT_DOUBLE_EQ(features.at("linked_share"), 1.0 / 3);

    // Without the genera the alternatives cannot be told apart: both shares are 0.
    profiler::MicrobialProfile without(*ref.loader);
    auto const& plain = fill(without);
    EXPECT_DOUBLE_EQ(plain.CongenerFitShare(), 0.0);
    EXPECT_DOUBLE_EQ(plain.OtherGenusFitShare(), 0.0);
    EXPECT_DOUBLE_EQ(plain.LowMapqShare(), 0.4);
    EXPECT_DOUBLE_EQ(plain.LinkedShare(), 1.0 / 3);
}

TEST(MicrobialProfile, DivergenceBeyondTheBaseQualities) {
    // ReadExcess: the differences per aligned base less the mean error probability of the record's bases.
    SamEntry sam;
    sam.m_cigar = "10M";
    sam.m_qual = "*";
    EXPECT_FALSE(profiler::ReadExcess(sam).has_value());
    sam.m_qual = std::string(10, 'I');  // Q40: 0.0001 per base
    EXPECT_NEAR(*profiler::ReadExcess(sam), -1e-4, 1e-7);
    sam.m_cigar = "8M2X";
    sam.m_qual = std::string(10, '+');  // Q10: 0.1
    EXPECT_NEAR(*profiler::ReadExcess(sam), 0.1, 1e-6);
    sam.m_cigar = "2S4M1D4M";  // a difference in 9 aligned bases; the clipped bases' qualities count too
    sam.m_qual = std::string(10, '5');  // Q20: 0.01
    EXPECT_NEAR(*profiler::ReadExcess(sam), 1.0 / 9 - 0.01, 1e-6);
    sam.m_cigar = "10S";
    EXPECT_FALSE(profiler::ReadExcess(sam).has_value());

    // A taxon's excess_median and excess_high_share, over all its best records (NoteRecord) with qualities.
    TinyReference ref;
    auto record = [&](std::string cigar, std::string qual) {
        SamEntry r;
        r.m_qname = "r";
        r.m_rname = "1_1";
        r.m_pos = 1;
        r.m_mapq = 60;
        r.m_cigar = std::move(cigar);
        r.m_seq = ref.gene.substr(0, 10);
        r.m_qual = std::move(qual);
        return r;
    };
    profiler::MicrobialProfile profile(*ref.loader);
    for (auto const& r : { record("8M2X", std::string(10, '+')),     // 0.1
                           record("10M", std::string(10, 'I')),      // -0.0001
                           record("9M1X", std::string(10, '5')) }) { // 0.09
        profile.NoteRecord(1, r);
        EXPECT_TRUE(profile.AddSam(1, 1, r, 1.0));
    }
    profile.NoteRecord(1, record("10M", "*"));  // no qualities: no excess
    profile.ApplyRecordEvidence();
    auto const& taxon = profile.GetTaxa().at(1);
    EXPECT_NEAR(taxon.ExcessMedian(), 0.09, 1e-6);
    EXPECT_DOUBLE_EQ(taxon.ExcessHighShare(), 2.0 / 3);
    std::map<std::string, double> features;
    for (auto const& [name, value] : profiler::TaxonFeatures(taxon)) features[name] = value;
    EXPECT_NEAR(features.at("excess_median"), 0.09, 1e-6);
    EXPECT_DOUBLE_EQ(features.at("excess_high_share"), 2.0 / 3);
    // Without the database's gene conservation factors the conservation pattern says nothing.
    EXPECT_DOUBLE_EQ(features.at("conserved_fast_depth_ratio"), 0.0);
    EXPECT_DOUBLE_EQ(features.at("conserved_hit_share"), 0.5);

    // Reads without qualities: both 0.
    profiler::MicrobialProfile unqualified(*ref.loader);
    auto plain = record("10M", "*");
    unqualified.NoteRecord(1, plain);
    plain.m_qual = std::string(10, 'I');
    EXPECT_TRUE(unqualified.AddSam(1, 1, plain, 1.0));
    unqualified.ApplyRecordEvidence();
    EXPECT_DOUBLE_EQ(unqualified.GetTaxa().at(1).ExcessMedian(), 0.0);
    EXPECT_DOUBLE_EQ(unqualified.GetTaxa().at(1).ExcessHighShare(), 0.0);
}

TEST(ProfileSam, ALongReadsGenesAreOneLinkedRead) {
    // A 200 bp gene, so that alignments pass the profiler's minimum length (more than 50 bases).
    std::string gene;
    uint32_t state = 7;
    for (int i = 0; i < 200; i++) {
        state = state * 1103515245u + 12345u;
        gene += "ACGT"[(state >> 16) & 3];
    }
    ScratchDir dir;
    std::string header = ">1_1\n";
    auto fna = dir.Write("reference.fna", header + gene + '\n');
    auto map = dir.Write("reference.map", "1\t1\t" + std::to_string(header.size()) + '\t' + std::to_string(header.size() + gene.size()) + '\n');
    GenomeLoader loader(fna, map);
    loader.LoadAllGenomes();

    std::string const tags = "\tZU:i:1\tZT:i:0\tZA:Z:*";
    // Read c: two parts on the gene (a primary and a supplementary record); read d: one part; read e: MAPQ 2,
    // which the profiler's MAPQ filter leaves out of the hits.
    auto sam = dir.Write("long.sam", "@HD\tVN:1.6\n" + kSamReadTypeComment + "pb\n" +
            Record("c", 0, "1_1", "60M", gene.substr(0, 60), tags, 1) +
            Record("c", 0x800, "1_1", "60M", gene.substr(100, 60), tags, 101) +
            Record("d", 0, "1_1", "60M", gene.substr(20, 60), tags, 21) +
            Record("e", 0, "1_1", "60M", gene.substr(40, 60), tags, 41, 2));
    profiler::Profiler profiler(loader);
    profiler::MicrobialProfile profile(loader);
    EXPECT_EQ(profiler.ProfileSam(sam, profile), "");
    ASSERT_EQ(profile.GetTaxa().size(), 1u);
    auto const& taxon = profile.GetTaxa().at(1);
    EXPECT_EQ(taxon.TotalHits(), 3u);
    EXPECT_DOUBLE_EQ(taxon.LinkedShare(), 0.5);
    EXPECT_DOUBLE_EQ(taxon.LowMapqShare(), 0.25);  // read e counts: 1 of the 4 records
}
