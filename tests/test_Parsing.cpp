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
#include "TestReference.h"

using namespace protal;
namespace fs = std::filesystem;
using protal::test::ScratchDir;

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
    struct MapLists {
        std::string output_dir, strain_dir, misc_dir;
        std::vector<std::string> prefixes, firsts, seconds, sams, profiles, names, truths, read_types, unmapped;

        bool Load(std::string const& path) {
            return Options::LoadFromMap(path, output_dir, strain_dir, misc_dir, prefixes, firsts, seconds, sams, profiles, names, truths, read_types, unmapped);
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

// A relative output folder (-o, else #OUTPUT_DIR) is relative to the folder protal runs in, and the folders in it are named
// once: up to 2026-10-08 a run wrote rel/rel/profiles (the folders the map did not name; named ones and the files named
// after the prefixes were right).
TEST(SampleMap, ARelativeOutputFolderIsUsedOnce) {
    ScratchDir dir;
    auto const map = dir.Write("rel.map", "#OUTPUT_DIR\tignored\n#SAM_OUTPUT_DIR\taln\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tSAM\tPROFILE\n"
                                          "a\ta\ta_1.fq\ta_2.fq\ta.sam.zst\ta.profile\n");
    auto const before = fs::current_path();
    fs::current_path(dir.path);  // the map's folders are made where protal runs
    MapLists lists;
    lists.output_dir = "rel";  // -o
    bool const ok = lists.Load(map);
    MapLists by_map;  // the map's #OUTPUT_DIR, without SAM and PROFILE columns
    bool const ok_by_map = by_map.Load(dir.Write("prefix.map", "#OUTPUT_DIR\tbymap\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\nb\tb\tb_1.fq\tb_2.fq\n"));
    fs::current_path(before);
    ASSERT_TRUE(ok);
    EXPECT_EQ(lists.sams, (Tokens{ "rel/aln/a.sam.zst" }));
    EXPECT_EQ(lists.profiles, (Tokens{ "rel/profiles/a.profile" }));
    EXPECT_EQ(lists.prefixes, (Tokens{ "rel/a" }));
    EXPECT_EQ((std::pair{ lists.strain_dir, lists.misc_dir }), (std::pair{ std::string("rel/strains"), std::string("rel/misc") }));
    EXPECT_TRUE(fs::is_directory(dir.path / "rel" / "profiles"));
    EXPECT_FALSE(fs::exists(dir.path / "rel" / "rel"));
    EXPECT_FALSE(fs::exists(dir.path / "ignored"));
    ASSERT_TRUE(ok_by_map);
    EXPECT_EQ(by_map.prefixes, (Tokens{ "bymap/b" }));
    EXPECT_EQ(by_map.sams, (Tokens{ "bymap/b.sam.zst" }));  // named after the prefix, in the output folder
    EXPECT_EQ(by_map.profiles, (Tokens{ "bymap/b.profile" }));
    EXPECT_EQ(by_map.misc_dir, "bymap/misc");
    EXPECT_FALSE(fs::exists(dir.path / "bymap" / "bymap"));
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

TEST(SampleMap, UnmappedReadsChoosePerSample) {
    ScratchDir dir;
    auto const header = "#OUTPUT_DIR\t" + (dir.path / "out").string() + "\n#SAMPLEID\tPREFIX\tFIRST\tSECOND\tUNMAPPED_READS\n";
    MapLists lists;
    ASSERT_TRUE(lists.Load(dir.Write("unmapped.map", header + "a\ta\ta_1.fq\ta_2.fq\twrite\nb\tb\tb_1.fq\tb_2.fq\tcount\n"
                                                            "c\tc\tc_1.fq\tc_2.fq\t-\n")));
    EXPECT_EQ(lists.unmapped, (Tokens{ "write", "count", "-" }));
    // '-' is --write_unmapped_reads's choice.
    EXPECT_EQ(Options::UnmappedReadsFlags(lists.unmapped, false), (std::vector<char>{ 1, 0, 0 }));
    EXPECT_EQ(Options::UnmappedReadsFlags(lists.unmapped, true), (std::vector<char>{ 1, 0, 1 }));
    MapLists none;
    ASSERT_TRUE(none.Load(dir.Write("none.map", MapHeader(dir.path / "out") + "a\ta\ta_1.fq\ta_2.fq\ta.profile\n")));
    EXPECT_TRUE(none.unmapped.empty());

    // Per sample through Options; --full_sam_header writes them for every sample.
    OptionsData d;
    d.unmapped_reads_list = lists.unmapped;
    Options by_map(d);
    EXPECT_TRUE(by_map.WriteUnmappedReads(0));
    EXPECT_FALSE(by_map.WriteUnmappedReads(1));
    EXPECT_FALSE(by_map.WriteUnmappedReads(2));
    EXPECT_FALSE(by_map.WriteUnmappedReads());
    d.write_unmapped_reads = true;
    Options all(d);
    EXPECT_FALSE(all.WriteUnmappedReads(1));
    EXPECT_TRUE(all.WriteUnmappedReads(2));
    d.full_sam_header = true;
    EXPECT_TRUE(Options(d).WriteUnmappedReads(1));

    testing::internal::CaptureStderr();
    MapLists bad;
    EXPECT_FALSE(bad.Load(dir.Write("bad.map", header + "a\ta\ta_1.fq\ta_2.fq\tyes\n")));
    auto const log = testing::internal::GetCapturedStderr();
    EXPECT_NE(log.find("UNMAPPED_READS is 'yes'"), std::string::npos) << log;
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
    struct TinyReference : test::LoadedReference {
        std::string const gene = test::TinyReference::kGene1;  // 50 bp

        TinyReference() : LoadedReference({ { 1, { test::TinyReference::kGene1 } } }, "parsing reference") {}
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

// Without a read type, a single read file's first reads decide (ReadTypeDetection.h) for the samples of the run, with
// a note; a read type given wins, and the reads of samples outside the run are not looked at.
TEST(Options, TheReadTypeOfASingleReadFileFromItsReads) {
    ScratchDir dir;
    std::string reads;
    for (int i = 0; i < 50; i++) reads += "@r" + std::to_string(i) + "\n" + std::string(4000, 'A') + "\n+\n" + std::string(4000, '3') + "\n";
    std::string const ont = dir.Write("ont.fq", reads);
    OptionsData data;
    data.first_list = { ont, ont, ont };
    data.second_list = { "", "", "" };
    data.prefix_list = { "a", "b", "c" };
    data.read_type_list = { "", "se", "" };
    data.range = { 0, 1 };
    Options options(data);
    std::vector<std::string> warnings, notes;
    options.ResolveReadTypes(warnings, &notes);
    EXPECT_EQ(options.GetReadType(0), ReadType::ONT);
    EXPECT_EQ(options.GetReadType(1), ReadType::Single);  // given
    EXPECT_EQ(options.GetReadType(2), ReadType::Single);  // not in the run
    ASSERT_EQ(notes.size(), 1u);
    EXPECT_NE(notes[0].find("Sample a: ONT reads (the first 50 reads up to 4.0 kb long, median 4.0 kb, median read quality Q18.0"),
              std::string::npos) << notes[0];
    EXPECT_NE(notes[0].find("aligned and profiled as such; --read_type (or a map's READ_TYPE) sets the kind of reads"), std::string::npos);
    EXPECT_TRUE(warnings.empty());
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
            profile.NoteRecord(1, 1, sam);
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
        profile.NoteRecord(1, 1, sam_at(11, 0, "2:0"));
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
    // A taxon's excess_median and excess_high_share, over all its best records (NoteRecord) with qualities: each record's
    // ReadExcess, its differences per aligned base less the mean error probability of its bases.
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
        profile.NoteRecord(1, 1, r);
        EXPECT_TRUE(profile.AddSam(1, 1, r, 1.0));
    }
    profile.NoteRecord(1, 1, record("10M", "*"));  // no qualities: no excess
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
    unqualified.NoteRecord(1, 1, plain);
    plain.m_qual = std::string(10, 'I');
    EXPECT_TRUE(unqualified.AddSam(1, 1, plain, 1.0));
    unqualified.ApplyRecordEvidence();
    EXPECT_DOUBLE_EQ(unqualified.GetTaxa().at(1).ExcessMedian(), 0.0);
    EXPECT_DOUBLE_EQ(unqualified.GetTaxa().at(1).ExcessHighShare(), 0.0);
}

TEST(MicrobialProfile, ALongReadsGenesAreOneLinkedRead) {
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

// NormalizeCigar against CIGARs normalised by hand. '=' becomes M, runs of one op merge (also where '=' makes them one),
// leading zeros go, and hard clips go, counted at the start and at the end. Rejected: an empty or malformed CIGAR, a zero
// count, an op the profiler does not walk (N, P, ...), nothing but hard clips, or not SEQ's length.
TEST(SamParsing, NormalizeCigarGivesTheNormalForm) {
    struct Case {
        std::string cigar;
        size_t length;
        char const* normal;  // nullptr: rejected
        uint32_t clip_start = 0, clip_end = 0;
    };
    std::vector<Case> const cases = {
            { "4M", 4, "4M" },
            { "2S3M1X2I1D4M3S", 15, "2S3M1X2I1D4M3S" },  // a deletion takes no base of SEQ
            { "2I2D2M", 4, "2I2D2M" },                    // an insertion and a deletion side by side stay apart
            { "2=1X1=", 4, "2M1X1M" },
            { "2=2M", 4, "4M" },
            { "1S2M1M", 4, "1S3M" },
            { "1X1X2M", 4, "2X2M" },
            { "2S1S2M", 5, "3S2M" },
            { "04M", 4, "4M" },
            { "000000000004M", 4, "4M" },
            { "1000000000M", 1000000000, "1000000000M" },  // more digits than the fast path reads
            { "5H2=1X1=3H", 4, "2M1X1M", 5, 3 },
            { "5H4M", 4, "4M", 5, 0 },
            { "4M3H", 4, "4M", 0, 3 },
            { "2H3H4M", 4, "4M", 5, 0 },
            { "2M3H2M", 4, "4M", 0, 0 },                   // between other ops: at neither end
            { "99999999999H4M", 4, "4M", UINT32_MAX, 0 },  // clips beyond 32 bits are capped
            { "", 0, nullptr },
            { "4M", 5, nullptr },
            { "2M1D2M", 5, nullptr },
            { "0M4M", 4, nullptr },
            { "4", 4, nullptr },
            { "M", 1, nullptr },
            { "*", 4, nullptr },
            { "4m", 4, nullptr },
            { "2M2P2M", 4, nullptr },
            { "2M5N2M", 4, nullptr },
            { "5H", 0, nullptr },
            { "99999999999999999999M", 4, nullptr },  // beyond 64 bits
    };
    for (auto const& c : cases) {
        std::string cigar = c.cigar;
        uint32_t clip_start = 7, clip_end = 7;
        bool const ok = NormalizeCigar(cigar, c.length, &clip_start, &clip_end);
        if (!c.normal) {
            EXPECT_FALSE(ok) << c.cigar;
            EXPECT_EQ(cigar, c.cigar) << "a rejected CIGAR is left as it was";
            continue;
        }
        ASSERT_TRUE(ok) << c.cigar;
        EXPECT_EQ(cigar, c.normal) << c.cigar;
        EXPECT_EQ(clip_start, c.clip_start) << c.cigar;
        EXPECT_EQ(clip_end, c.clip_end) << c.cigar;
    }
}

namespace {
    struct Lcg {
        uint64_t state;
        uint32_t Next(uint32_t n) {
            state = state * 6364136223846793005ull + 1442695040888963407ull;
            return static_cast<uint32_t>((state >> 33) % n);
        }
    };

    // A CIGAR's ops one per column, '=' as M, without its hard clips; empty if it has more than max_columns.
    std::string Columns(std::string const& cigar, uint64_t max_columns = 10000) {
        std::string columns;
        uint64_t count = 0;
        for (char const c : cigar) {
            if (c >= '0' && c <= '9') {
                count = count * 10 + static_cast<uint64_t>(c - '0');
                if (count > max_columns) return {};
                continue;
            }
            if (c != 'H') columns.append(count, c == '=' ? 'M' : c);
            if (columns.size() > max_columns) return {};
            count = 0;
        }
        return columns;
    }
}

// On random CIGARs, mostly well-formed: one IsNormalCigar takes is left as it is (the fast path); one NormalizeCigar
// rewrites has the same columns but for hard clips, and is left as it is when normalised again.
TEST(SamParsing, NormalizeCigarKeepsTheColumnsAndANormalCigar) {
    Lcg random{ 17 };
    std::string const ops = "MXIDS=HNP*";
    size_t fast = 0, rewritten = 0;
    for (int n = 0; n < 50000; n++) {
        std::string cigar;
        size_t query = 0;
        int const parts = 1 + static_cast<int>(random.Next(6));
        for (int p = 0; p < parts; p++) {
            uint32_t const kind = random.Next(20);
            std::string count = std::to_string(kind == 0 ? 0 : 1 + random.Next(kind < 3 ? 2000000000u : 150));
            if (kind == 1) count = "0" + count;           // a leading zero
            if (kind == 2) count = std::string(12, '9');  // longer than the fast path reads
            char const op = random.Next(10) < 7 ? "MXIDS"[random.Next(5)] : ops[random.Next(static_cast<uint32_t>(ops.size()))];
            if (random.Next(30) != 0) cigar += count;
            if (random.Next(40) != 0) cigar += op;
            if (op != 'D' && op != 'H') query += std::stoul(count);
        }
        // The read's length right, or off.
        size_t const length = random.Next(4) == 0 ? query + random.Next(3) : query;
        std::string normal = cigar;
        uint32_t clip_start = 7, clip_end = 7;
        bool const ok = NormalizeCigar(normal, length, &clip_start, &clip_end);
        if (IsNormalCigar(cigar, length)) {
            fast++;
            ASSERT_TRUE(ok) << cigar;
            ASSERT_EQ(normal, cigar);
            ASSERT_EQ(clip_start, 0u) << cigar;
            ASSERT_EQ(clip_end, 0u) << cigar;
            continue;
        }
        if (!ok) {
            ASSERT_EQ(normal, cigar) << "a rejected CIGAR is left as it was";
            continue;
        }
        rewritten++;
        ASSERT_EQ(normal.find_first_of("=HNP*"), std::string::npos) << cigar << " -> " << normal;
        if (auto const columns = Columns(cigar); !columns.empty()) ASSERT_EQ(Columns(normal), columns) << cigar;
        std::string again = normal;
        ASSERT_TRUE(NormalizeCigar(again, length, &clip_start, &clip_end)) << normal;
        ASSERT_EQ(again, normal);
        ASSERT_EQ(clip_start, 0u) << normal;
        ASSERT_EQ(clip_end, 0u) << normal;
    }
    EXPECT_GT(fast, 1000u);
    EXPECT_GT(rewritten, 1000u);
}

TEST(SamParsing, FieldsAreSplitAsLineSplitterSplitsThem) {
    Lcg random{ 5 };
    std::vector<std::string> tokens;
    std::vector<std::string_view> fields;
    for (int n = 0; n < 20000; n++) {
        std::string line;
        int const length = static_cast<int>(random.Next(30));
        for (int i = 0; i < length; i++) line += "ab\t\t1_2*"[random.Next(8)];
        LineSplitter::Split(line, "\t", tokens);
        sam_detail::SplitFields(line, fields);
        ASSERT_EQ(fields.size(), tokens.size()) << line;
        for (size_t i = 0; i < fields.size(); i++) ASSERT_EQ(std::string(fields[i]), tokens[i]) << line;
    }
}

TEST(SamParsing, CountsAndNameNumbersParseAsStoiAndStoul) {
    for (std::string const s : { "0", "7", "150", "2147483647", "000123", "4294967295", "18446744073709551615" }) {
        EXPECT_EQ(NameNumber(s, 0, s.size()), std::stoul(s)) << s;
    }
    std::string const name = "123_4567_rest";
    EXPECT_EQ(NameNumber(name, 0, 3), 123u);
    EXPECT_EQ(NameNumber(name, 4, 4), 4567u);
    std::string const cigar = "12M3X2147483647S";
    EXPECT_EQ(CigarCount(cigar, 0, 2), 12);
    EXPECT_EQ(CigarCount(cigar, 3, 4), 3);
    EXPECT_EQ(CigarCount(cigar, 5, 15), 2147483647);  // the largest int
    CigarInfo info;
    CompressedCigarInfo("5S20M2I3D10M1X", info);
    EXPECT_EQ(info.softclipped, 5);
    EXPECT_EQ(info.matches, 30);
    EXPECT_EQ(info.insertions, 2);
    EXPECT_EQ(info.deletions, 3);
    EXPECT_EQ(info.mismatches, 1);
    EXPECT_EQ(info.clipped_alignment_length, 36);
    EXPECT_EQ(AlignmentLengthRef("5S20M2I3D10M1X"), 34u);
}

TEST(SamReader, ATextReadsAsTheSameStream) {
    std::string const text = "@HD\tVN:1.6\r\n\n" + Record("a", kPaired | kBothAlign | kRead1, "1_1", "4M", "ACGT") +
                             Record("a", kPaired | kBothAlign | kRead2, "1_1", "2M1X1M", "ACGA", "\tZU:i:3\tZA:Z:2:1") +
                             "u\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\r\n" + Record("b", 16, "1_2", "3=1M", "ACGT", "") + "c\t0\t1_1";
    std::istringstream in(text);
    std::vector<std::string> headers_stream, headers_text;
    SamReader stream(in, [&](std::string const& line) { headers_stream.push_back(line); });
    SamReader view(std::string_view(text), [&](std::string const& line) { headers_text.push_back(line); }, 0);
    SamEntry s1, s2, t1, t2;
    bool hs1 = false, hs2 = false, ht1 = false, ht2 = false;
    std::string stream_error, text_error;
    while (true) {
        bool more_stream = false, more_text = false;
        try { more_stream = stream.Next(s1, s2, hs1, hs2); } catch (SamFormatError const& e) { stream_error = e.what(); }
        try { more_text = view.Next(t1, t2, ht1, ht2); } catch (SamFormatError const& e) { text_error = e.what(); }
        ASSERT_EQ(more_stream, more_text);
        if (!more_stream) break;
        ASSERT_EQ(hs1, ht1);
        ASSERT_EQ(hs2, ht2);
        if (hs1) EXPECT_EQ(s1.ToString(), t1.ToString());
        if (hs2) EXPECT_EQ(s2.ToString(), t2.ToString());
    }
    EXPECT_EQ(stream_error, text_error);
    EXPECT_NE(stream_error.find("line 7"), std::string::npos) << stream_error;
    EXPECT_EQ(headers_stream, headers_text);
    EXPECT_EQ(headers_text, std::vector<std::string>{ "@HD\tVN:1.6" });
    EXPECT_EQ(stream.Records(), view.Records());
    EXPECT_EQ(stream.Skipped(), view.Skipped());
    EXPECT_EQ(stream.RecordsWithoutTags(), view.RecordsWithoutTags());
    // A text that is a part of a file numbers its lines from there.
    std::string const bad = "x\t0\n";
    SamReader later(std::string_view(bad), {}, 41);
    SamEntry e1, e2;
    bool h1 = false, h2 = false;
    try {
        later.Next(e1, e2, h1, h2);
        ADD_FAILURE() << "no error";
    } catch (SamFormatError const& e) {
        EXPECT_NE(std::string(e.what()).find("line 42: expected at least 11 tab-separated fields, found 2"), std::string::npos) << e.what();
    }
}

namespace {
    // What loading the sample map at `path` says, failing the test if the map loads.
    std::string MapError(std::string const& path) {
        testing::internal::CaptureStderr();
        MapLists lists;
        bool const loaded = lists.Load(path);
        auto const log = testing::internal::GetCapturedStderr();
        EXPECT_FALSE(loaded) << path;
        return log;
    }
}

TEST(SampleMap, AMissingFileIsAnError) {
    ScratchDir dir;
    auto const log = MapError((dir.path / "none.map").string());
    EXPECT_NE(log.find("none.map does not exist"), std::string::npos) << log;
}

TEST(SampleMap, AMalformedHeaderIsAnError) {
    ScratchDir dir;
    auto const out = "#OUTPUT_DIR\t" + (dir.path / "out").string() + "\n";
    auto expect = [&](std::string const& name, std::string const& map, std::string const& message) {
        auto const log = MapError(dir.Write(name, map));
        EXPECT_NE(log.find(message), std::string::npos) << name << ": " << log;
    };
    expect("twice.map", out + "#SAMPLEID\tPREFIX\tFIRST\tFIRST\ns\ts\ta.fq\tb.fq\n", "Column 'FIRST' is defined twice");
    // A '#' line among the rows, after the column header.
    expect("late.map", out + "#SAMPLEID\tPREFIX\tFIRST\ns\ts\ta.fq\n# a comment\n",
           "Line 4: Did not expect header line but line starts with #");
    expect("no_output.map", "#SAMPLEID\tPREFIX\tFIRST\ns\ts\ta.fq\n", "Line 1: Output directory not defined");
    expect("no_value.map", "#OUTPUT_DIR\n#SAMPLEID\tPREFIX\tFIRST\ns\ts\ta.fq\n", "Line 1: Expected value for key #OUTPUT_DIR");
    // Rows before the column header, and a header without FIRST.
    expect("no_header.map", out + "s\ts\ta.fq\n", "The columns must be specified: PREFIX, FIRST");
    expect("no_first.map", out + "#SAMPLEID\tPREFIX\ns\ts\n", "The columns must be specified: PREFIX, FIRST");
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

// Sample IDs name files (misc/<sample>_runtime.tsv): what Linux and macOS cannot hold safely is refused.
TEST(Options, SampleIdsThatCannotNameFiles) {
    using Ids = std::vector<std::string>;
    for (auto const& id : Ids{ "S1", "sample_1.2-x", "Probe_\xc3\x84", std::string(Options::kMaxSampleIdBytes, 'x') }) {
        EXPECT_EQ(Options::SampleIdProblem(id), "") << id;
    }
    for (auto const& id : Ids{ "", ".", "..", "a/b", "a:b", "-a", "a\tb", "a\x01" "b", "a\xff" "b", "a\xc3",
                               std::string(Options::kMaxSampleIdBytes + 1, 'x') }) {
        EXPECT_NE(Options::SampleIdProblem(id), "") << id;
    }
}

// --profile_only's items: files as they are, wildcards expanded to the SAM files they match (sorted), errors for a
// pattern without SAMs and for a folder; a missing file is left for the existence check.
TEST(Options, SamPatternsExpandToTheSamFilesTheyMatch) {
    ScratchDir dir("sam patterns");
    for (auto const* folder : { "r2/alignments", "r1/alignments" }) fs::create_directories(dir.path / folder);
    auto const b = dir.Write("r1/alignments/b.sam.zst", "");
    auto const a = dir.Write("r1/alignments/a.sam", "");
    auto const c = dir.Write("r2/alignments/c.sam.gz", "");
    dir.Write("r1/alignments/a.sam.err", "");
    dir.Write("r1/alignments/b.sam.zst.partial", "");
    dir.Write("r1/alignments/a.profile", "");
    std::vector<std::string> errors, notes;
    auto const sams = Options::ExpandSamFiles({ dir / "r*/alignments/*", dir / "x.sam.zst" }, errors, notes);
    EXPECT_EQ(sams, (std::vector<std::string>{ a, b, c, dir / "x.sam.zst" }));
    EXPECT_TRUE(errors.empty());
    ASSERT_EQ(notes.size(), 1u);
    EXPECT_NE(notes[0].find("matches 3 SAM file(s) (and 3 other file(s) or folder(s), passed over)"), std::string::npos) << notes[0];

    // A file whose name holds a wildcard is that file.
    auto const odd = dir.Write("r1/alignments/odd[1].sam", "");
    EXPECT_EQ(Options::ExpandSamFiles({ odd }, errors, notes), (std::vector<std::string>{ odd }));

    errors.clear();
    EXPECT_TRUE(Options::ExpandSamFiles({ dir / "r*/alignments/*.err", dir / "r1", dir / "none/*.sam" }, errors, notes).empty());
    ASSERT_EQ(errors.size(), 3u);
    EXPECT_NE(errors[0].find("matches no SAM file"), std::string::npos) << errors[0];
    EXPECT_NE(errors[1].find("is a folder"), std::string::npos) << errors[1];
    EXPECT_NE(errors[2].find("matches no SAM file"), std::string::npos) << errors[2];
}

// The samples' names: the SAMs' file names, or, where they repeat, the folders in which the paths differ.
TEST(Options, SamSampleNamesFromFilesOrFolders) {
    using Names = std::vector<std::string>;
    EXPECT_EQ(Options::SamSampleNames({ "/r1/a.sam.zst", "/r2/b.sam.gz", "c.sam" }), (Names{ "a", "b", "c" }));
    EXPECT_EQ(Options::SamSampleNames({ "/d/S1/aln.sam.zst", "/d/S2/aln.sam.zst" }), (Names{ "S1", "S2" }));
    EXPECT_EQ(Options::SamSampleNames({ "/p/study1/alignments/sa.sam.zst", "/p/study2/alignments/sa.sam.zst",
                                        "/p/study2/alignments/sb.sam.zst" }),
              (Names{ "study1_sa", "study2_sa", "study2_sb" }));
    // Paths of different depths: every level from the shortest path's end differs.
    EXPECT_EQ(Options::SamSampleNames({ "/d/a/x.sam", "/d/b/c/x.sam" }), (Names{ "a", "b_c" }));
    // The same file twice keeps one name, for the duplicate check to report.
    EXPECT_EQ(Options::SamSampleNames({ "/d/a/x.sam", "/d/a/x.sam" }), (Names{ "x", "x" }));
}
