// The reference genomes of simulate_metagenomes's read simulators (src/RandomForest/GenomeStore.h): files read whole
// and inflated (gzip of one or several members, BGZF, zstd, plain), FASTA records parsed in bulk exactly as the
// line-by-line reader before them took them, the genome store's 2-bit files (any slice the same as the sequence's,
// other characters kept, a changed FASTA read again), contigs loaded alike with and without the store, and a genome
// table whose first row is data although its fields hold a header's words.

#include <gtest/gtest.h>

#include <cctype>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <zlib-ng.h>
#include <zstd.h>

#include "IO/Bgzf.h"
#include "RandomForest/GenomeStore.h"
#include "RandomForest/MetagenomeSimulator.h"
#include "RandomForest/ReadPipeline.h"
#include "TestUtil.h"

using namespace protal::sim;
using protal::test::ScratchDir;
namespace fs = std::filesystem;

namespace {

// The records as LoadContigs read them before the bulk parser (getline, then character by character).
FastaRecords LineByLine(std::string const& text) {
    FastaRecords records;
    std::istringstream in(text);
    std::string line;
    bool in_record = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            auto const end = line.find_first_of(" \t\r", 1);
            records.names.push_back(line.substr(1, end == std::string::npos ? std::string::npos : end - 1));
            records.seqs.emplace_back();
            in_record = true;
            continue;
        }
        if (!in_record) continue;
        for (char c : line) {
            if (!std::isspace(static_cast<unsigned char>(c))) {
                records.seqs.back() += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
            }
        }
    }
    return records;
}

void ExpectSameRecords(FastaRecords const& got, FastaRecords const& want, std::string const& what) {
    EXPECT_EQ(got.names, want.names) << what;
    ASSERT_EQ(got.seqs.size(), want.seqs.size()) << what;
    for (std::size_t i = 0; i < want.seqs.size(); ++i) EXPECT_EQ(got.seqs[i], want.seqs[i]) << what << " record " << i;
}

// One gzip member (zlib-ng); with full_header every optional header field (extra, name, comment, header CRC).
std::string GzipBytes(std::string_view text, bool full_header = false) {
    zng_stream zs{};
    EXPECT_EQ(zng_deflateInit2(&zs, 6, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY), Z_OK);
    static unsigned char extra[] = {'X', 'Y', 4, 0, 1, 2, 3, 4};
    static unsigned char name[] = "genome.fna", comment[] = "a comment";
    zng_gz_header header{};
    header.extra = extra;
    header.extra_len = sizeof(extra);
    header.name = name;
    header.comment = comment;
    header.hcrc = 1;
    if (full_header) EXPECT_EQ(zng_deflateSetHeader(&zs, &header), Z_OK);
    std::string out(zng_deflateBound(&zs, text.size()) + 64, '\0');
    zs.next_in = reinterpret_cast<uint8_t const*>(text.data());
    zs.avail_in = static_cast<uint32_t>(text.size());
    zs.next_out = reinterpret_cast<uint8_t*>(out.data());
    zs.avail_out = static_cast<uint32_t>(out.size());
    EXPECT_EQ(zng_deflate(&zs, Z_FINISH), Z_STREAM_END);
    out.resize(zs.total_out);
    zng_deflateEnd(&zs);
    return out;
}

std::string ZstdBytes(std::string_view text) {
    std::string out(ZSTD_compressBound(text.size()), '\0');
    std::size_t const n = ZSTD_compress(out.data(), out.size(), text.data(), text.size(), 3);
    EXPECT_FALSE(ZSTD_isError(n));
    out.resize(n);
    return out;
}

std::string BgzfBytes(std::string const& text) {
    std::string out;
    EXPECT_TRUE(protal::bgzf::Compress(text.data(), text.size(), out));
    out.append(reinterpret_cast<char const*>(protal::bgzf::kEof), sizeof(protal::bgzf::kEof));
    return out;
}

// A FASTA with the quirks genome files have: text before the first header, descriptions, CRLF lines, lower case,
// white space inside lines, empty lines and records, IUPAC codes, bytes above 127, no newline at the end.
std::string QuirkyFasta(std::mt19937& gen) {
    std::string text = "a line before any header\n>c1 a description\nacgtNNnn\r\nAC GT\tac\n\n>c2\tx\n>c3\r\nRYKMswbdhv\n";
    text += ">c4\nAC\xC3\xA9GT*-.\n>big\n";
    std::string const bases = "ACGTacgtNn";
    for (int line = 0; line < 400; ++line) {
        for (int i = 0; i < 80; ++i) text += bases[gen() % (line % 50 == 0 ? bases.size() : 4)];
        if (line % 97 == 0) text += " \t";
        text += line % 13 == 0 ? "\r\n" : "\n";
    }
    text += ">last no newline\nACGTTGCA";
    return text;
}

}  // namespace

TEST(GenomeStore, FastaParsedAsLineByLine) {
    std::mt19937 gen(3);
    std::string const quirky = QuirkyFasta(gen);
    ExpectSameRecords(ParseFasta(quirky), LineByLine(quirky), "quirky");
    for (std::string const text : {"", ">", ">\n", "no header at all\nACGT\n", ">a\n>b\n\n", ">x y\r\n\r\n acgt \r\n",
                                   "\n\n>only\nNNNN"}) {
        ExpectSameRecords(ParseFasta(text), LineByLine(text), "'" + text + "'");
    }
    // a genome-sized record of 80-base lines, as in the build's FASTAs
    std::string genome = ">contig1 desc\n";
    std::string const seq = protal::test::RandomSequence(200000, gen);
    for (std::size_t at = 0; at < seq.size(); at += 80) genome += seq.substr(at, 80) + "\n";
    auto const records = ParseFasta(genome);
    ASSERT_EQ(records.seqs.size(), 1u);
    EXPECT_EQ(records.names[0], "contig1");
    EXPECT_EQ(records.seqs[0], seq);
}

TEST(GenomeStore, WholeFilesInflated) {
    ScratchDir dir("genome store");
    std::mt19937 gen(5);
    std::string const text = QuirkyFasta(gen) + "\n>more\n" + protal::test::TestData(300000);
    EXPECT_EQ(ReadWholeFile(dir.Write("plain.fna", text)), text);
    EXPECT_EQ(ReadWholeFile(dir.Write("one.fna.gz", GzipBytes(text))), text);
    EXPECT_EQ(ReadWholeFile(dir.Write("header.fna.gz", GzipBytes(text, true))), text);
    std::string const half = text.substr(0, text.size() / 2), rest = text.substr(text.size() / 2);
    EXPECT_EQ(ReadWholeFile(dir.Write("two.fna.gz", GzipBytes(half) + GzipBytes(rest))), text);
    EXPECT_EQ(ReadWholeFile(dir.Write("padded.fna.gz", GzipBytes(text) + std::string(512, '\0'))), text);
    EXPECT_EQ(ReadWholeFile(dir.Write("bgzf.fna.gz", BgzfBytes(text))), text);
    EXPECT_EQ(ReadWholeFile(dir.Write("frames.fna.zst", ZstdBytes(half) + ZstdBytes(rest))), text);
    EXPECT_EQ(ReadWholeFile(dir.Write("named_plain.fna", ZstdBytes(text))), text);  // told by its bytes, not its name
    EXPECT_EQ(ReadWholeFile(dir.Write("empty.fna", "")), "");

    std::string const gz = GzipBytes(text), zst = ZstdBytes(text);
    for (auto const& [name, bytes] : std::vector<std::pair<std::string, std::string>>{
             {"cut.fna.gz", gz.substr(0, gz.size() - 9)}, {"cut_header.fna.gz", gz.substr(0, 6)},
             {"garbage.fna.gz", gz + "trailing text"}, {"cut.fna.zst", zst.substr(0, zst.size() - 5)}}) {
        EXPECT_THROW(ReadWholeFile(dir.Write(name, bytes)), std::runtime_error) << name;
    }
    EXPECT_THROW(ReadWholeFile(dir / "missing.fna"), std::runtime_error);
    std::string corrupt = gz;
    corrupt[corrupt.size() - 6] ^= 0x55;  // the trailer's CRC
    EXPECT_THROW(ReadWholeFile(dir.Write("crc.fna.gz", corrupt)), std::runtime_error);
}

TEST(GenomeStore, FilesGiveEverySlice) {
    ScratchDir dir("genome store");
    std::mt19937 gen(7);
    FastaRecords records;
    records.names = {"c1", "c2", "empty", "c4"};
    std::string c1 = protal::test::RandomSequence(10001, gen);
    for (std::size_t at : {0u, 1u, 2u, 3u, 4u, 5000u, 10000u}) c1[at] = 'N';
    c1.replace(100, 50, std::string(50, 'N'));
    c1.replace(150, 7, "RYKMSWB");
    std::string const c2 = "NNNNACGTNNNN*-ACGT";
    records.seqs = {c1, c2, "", protal::test::RandomSequence(3, gen)};
    std::string const fasta = dir.Write("g.fna", ">x\nA\n");
    std::string error;
    ASSERT_TRUE(GenomeFile::Write(dir.path / "store", fasta, records, error)) << error;
    auto const file = GenomeFile::Open(dir.path / "store", fasta);
    ASSERT_TRUE(file);
    ASSERT_EQ(file->Contigs(), 4u);
    std::string slice;
    for (std::size_t k = 0; k < 4; ++k) {
        EXPECT_EQ(file->Name(k), records.names[k]);
        ASSERT_EQ(file->Length(k), records.seqs[k].size());
        std::string const& seq = records.seqs[k];
        for (int i = 0; i < 2000; ++i) {
            std::uint64_t const start = seq.empty() ? 0 : gen() % (seq.size() + 1);
            std::uint64_t const length = gen() % 400;
            file->Extract(k, start, length, slice);
            ASSERT_EQ(slice, seq.substr(start, length)) << k << " " << start << " " << length;
        }
        file->Extract(k, 0, seq.size() + 10, slice);
        EXPECT_EQ(slice, seq);
    }
    EXPECT_THROW(file->Extract(1, c2.size() + 1, 1, slice), std::out_of_range);

    // The FASTA changed (its size and time): the store's file is not taken; another FASTA has another file.
    dir.Write("g.fna", ">x\nAC\n");
    fs::last_write_time(fasta, fs::last_write_time(fasta) + std::chrono::seconds(5));
    EXPECT_FALSE(GenomeFile::Open(dir.path / "store", fasta));
    EXPECT_NE(GenomeFile::PathOf(dir.path / "store", fasta), GenomeFile::PathOf(dir.path / "store", dir / "h.fna"));
    // A cut file is not taken either.
    ASSERT_TRUE(GenomeFile::Write(dir.path / "store", fasta, records, error)) << error;
    ASSERT_TRUE(GenomeFile::Open(dir.path / "store", fasta));
    fs::resize_file(GenomeFile::PathOf(dir.path / "store", fasta), 200);
    EXPECT_FALSE(GenomeFile::Open(dir.path / "store", fasta));
}

TEST(GenomeStore, ContigsAlikeWithAndWithoutTheStore) {
    ScratchDir dir("genome store");
    std::mt19937 gen(11);
    std::string text;
    std::vector<std::string> seqs;
    for (int c = 0; c < 12; ++c) {
        seqs.push_back(protal::test::RandomSequence(c == 3 ? 90 : 500 + gen() % 5000, gen));
        if (c == 5) seqs.back().replace(40, 30, std::string(30, 'N'));
        text += ">contig" + std::to_string(c) + " d\n";
        for (std::size_t at = 0; at < seqs.back().size(); at += 70) text += seqs.back().substr(at, 70) + "\n";
    }
    std::string const fasta = dir.Write("genome.fna.gz", GzipBytes(text));
    auto const plain = LoadContigs("G", fasta, 150);
    auto const written = LoadContigs("G", fasta, 150, dir.path / "store");  // reads the FASTA, writes the store's file
    ASSERT_TRUE(fs::exists(GenomeFile::PathOf(dir.path / "store", fasta)));
    auto const mapped = LoadContigs("G", fasta, 150, dir.path / "store");   // from the store
    EXPECT_FALSE(plain->file);
    EXPECT_FALSE(written->file);
    EXPECT_TRUE(mapped->file);
    for (auto const* contigs : {written.get(), mapped.get()}) {
        EXPECT_EQ(contigs->names, plain->names);
        EXPECT_EQ(contigs->lengths, plain->lengths);
        EXPECT_EQ(contigs->starts, plain->starts);
        EXPECT_EQ(contigs->longest, plain->longest);
        std::string a, b;
        for (std::size_t k = 0; k < plain->names.size(); ++k) {
            for (int i = 0; i < 200; ++i) {
                std::uint64_t const start = gen() % (plain->lengths[k] + 1), length = gen() % 600;
                plain->Extract(k, start, length, a);
                contigs->Extract(k, start, length, b);
                ASSERT_EQ(a, b);
            }
        }
    }
    EXPECT_EQ(plain->names.size(), 11u);  // contig3 (90 bases) is shorter than 150
    EXPECT_EQ(plain->names[0], "contig0");
    EXPECT_THROW(LoadContigs("G", fasta, 100000, dir.path / "store"), std::runtime_error);
    EXPECT_THROW(LoadContigs("G", dir / "missing.fna", 150, dir.path / "store"), std::runtime_error);
}

TEST(GenomeStore, GenomeTableFirstRowIsDataThoughItHoldsHeaderWords) {
    ScratchDir dir("genome table");
    // A path with "taxonomy", "genomes" and "files" in it, as GTDB's release folders and a project folder may give.
    std::string const folder = (dir.path / "my_taxonomy" / "genomic_files_all" / "gtdb_genomes_all_r226").string();
    std::string const rows = "G1\td__Bacteria;p__P;c__C;o__O;f__F;g__Taxeobacter;s__Taxeobacter x\t" + folder +
                             "/G1.fna.gz\t1000\nG2\td__Bacteria;p__P;c__C;o__O;f__F;g__G;s__G y\t" + folder + "/G2.fna.gz\t2000\n";
    auto const genomes = read_genome_table(dir.Write("genomes.tsv", rows));
    ASSERT_EQ(genomes.size(), 2u);
    EXPECT_EQ(genomes[0].name, "G1");
    EXPECT_EQ(genomes[0].genome_length, 1000u);
    EXPECT_EQ(genomes[1].fasta_path, folder + "/G2.fna.gz");
    // a real header is still one
    auto const headed = read_genome_table(dir.Write("headed.tsv", "genome\ttaxonomy\tfasta_path\tgenome_length\n" + rows));
    ASSERT_EQ(headed.size(), 2u);
    EXPECT_EQ(headed[0].name, "G1");
    EXPECT_EQ(headed[1].genome_length, 2000u);
}
