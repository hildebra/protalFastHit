// Unit tests for parsing tab-separated input: the line splitter, sample map rows, and SAM records
// as the profiler reads them (protal's own and other aligners').
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
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
                                                  "p\tp\tp_1.fq\tp_2.fq\tshort\nl\tl\tl.fq.gz\t-\tpacbio\n")));
    EXPECT_EQ(typed.read_types, (Tokens{ "short", "pacbio" }));
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
        (type == ReadType::PacBio ? d.model_pacbio : d.model_se) = own;
        return fs::path(Options(d).ModelDbFile(type).Path()).filename().string();
    };
    EXPECT_EQ(model("", "", ReadType::Paired), "model.xml");
    EXPECT_EQ(model("", "", ReadType::Single), "model_se.xml");
    EXPECT_EQ(model("", "", ReadType::PacBio), "model_pacbio.xml");
    // --model replaces all, unless --model_se or --model_pacbio is given for their reads.
    EXPECT_EQ(model("other", "", ReadType::Paired), "other.xml");
    EXPECT_EQ(model("other", "", ReadType::Single), "other.xml");
    EXPECT_EQ(model("other", "", ReadType::PacBio), "other.xml");
    EXPECT_EQ(model("other", "se2.xml", ReadType::Single), "se2.xml");
    EXPECT_EQ(model("other", "se2.xml", ReadType::Paired), "other.xml");
    EXPECT_EQ(model("other", "hifi", ReadType::PacBio), "hifi.xml");
    auto const existing = dir.Write("elsewhere.xml", "");
    EXPECT_EQ(model("", existing, ReadType::Single), "elsewhere.xml");
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

TEST(Options, ReadTypesOfSamples) {
    ScratchDir dir;
    // From the read files: no second file means single-end reads.
    OptionsData reads;
    // PacBio reads come as READ_TYPE (or --read_type) pacbio.
    reads.first_list = { "p_1.fq", "s.fq", "l.fq" };
    reads.second_list = { "p_2.fq", "", "" };
    reads.prefix_list = { "p", "s", "l" };
    reads.read_type_list = { READ_TYPE_SHORT, READ_TYPE_SHORT, READ_TYPE_PACBIO };
    reads.range = { 0, 1 };
    Options from_reads(reads);
    from_reads.ResolveReadTypes();
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
                      dir.Write("long.sam", "@HD\tVN:1.6\n" + kSamReadsComment + "PacBio\n" + Record("c", 0x800, "1_1", "5H20M", seq)) };
    sams.prefix_list = { "pe", "se", "empty", "long" };
    sams.range = { 1 };
    Options from_sams(sams);
    from_sams.ResolveReadTypes();
    EXPECT_EQ(from_sams.GetReadType(0), ReadType::Paired);
    EXPECT_EQ(from_sams.GetReadType(1), ReadType::Single);
    EXPECT_EQ(from_sams.GetReadType(2), ReadType::Paired);
    EXPECT_EQ(from_sams.GetReadType(3), ReadType::PacBio);
    EXPECT_TRUE(from_sams.AnySample(ReadType::Single));
    EXPECT_FALSE(from_sams.AnySample(ReadType::Paired));  // only the single-end sample is in the range
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
