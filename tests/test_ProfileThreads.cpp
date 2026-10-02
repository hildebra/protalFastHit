// Profiling a SAM on several threads (Profiler::ProfileSam, SamChunks.h): the chunks a SAM is cut into, and
// a profile that is the same as on one thread, taxon for taxon and gene for gene, in the order of the maps.
#include <gtest/gtest.h>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include <unistd.h>
#include "Profiling/Profiler.h"
#include "Profiling/SamChunks.h"
#include "IO/SamFile.h"

using namespace protal;
namespace fs = std::filesystem;

namespace {
    struct ScratchDir {
        fs::path path;
        ScratchDir() {
            path = fs::temp_directory_path() / ("protal profile threads test " + std::to_string(::getpid()));
            fs::create_directories(path);
        }
        ~ScratchDir() { fs::remove_all(path); }

        std::string Write(std::string const& name, std::string const& content) const {
            auto file = path / name;
            std::ofstream(file, std::ios::binary) << content;
            return file.string();
        }
    };

    struct Random {
        uint64_t state;
        explicit Random(uint64_t seed) : state(seed) {}
        uint32_t Next(uint32_t n) {
            state = state * 6364136223846793005ull + 1442695040888963407ull;
            return static_cast<uint32_t>((state >> 33) % n);
        }
    };

    constexpr int kPaired = 0x1, kBothAlign = 0x2, kReverse = 0x10, kMateReverse = 0x20, kRead1 = 0x40, kRead2 = 0x80,
                  kSecondary = 0x100, kSupplementary = 0x800;
    constexpr int kTaxa = 5, kGenes = 4, kGeneLength = 300;

    // Gene neighbours of the taxa (families 100-102, order 200): in each clade gene g's 3' end faces gene g + 1's
    // 5' end (gene 4's, gene 1's) in some of the species informative there and nothing in the others, so that the
    // links of the long reads of Sam (their parts on genes g, g + 1, ...) have smoothed shares of many values.
    gene_neighbours::Table Neighbours() {
        std::string text = "# protal gene neighbours: genomes=90 species=60 max_gap=3000 ranks=family,order\n"
                           "clade\tgene\tend\tpartner\tpartner_end\tspecies\tinformative\tgap_median\tgap_min\tgap_max\n";
        for (int clade : { 100, 101, 102, 200 }) {
            for (int gene = 1; gene <= kGenes; gene++) {
                int const informative = 7 + (clade * 3 + gene * 5) % 11;
                int const species = 1 + (clade + gene * 7) % (informative - 1);
                std::string const end = std::to_string(clade) + '\t' + std::to_string(gene) + "\t3\t";
                std::string const of = '\t' + std::to_string(informative);
                text += end + std::to_string(gene % kGenes + 1) + "\t5\t" + std::to_string(species) + of + "\t20\t10\t30\n";
                text += end + "0\t0\t" + std::to_string(informative - species) + of + "\t0\t0\t0\n";
            }
        }
        gene_neighbours::Table table;
        std::istringstream is(text);
        EXPECT_EQ(table.Read(is), "");
        for (uint32_t taxid = 1; taxid <= kTaxa; taxid++) table.SetLineage(taxid, { taxid, 100 + (taxid - 1) / 2, 200, 300 });
        return table;
    }

    // kTaxa taxa (1..kTaxa) of kGenes genes (1..kGenes) each, of random bases; with_neighbours: with Neighbours().
    struct Reference {
        ScratchDir dir;
        std::map<std::pair<int, int>, std::string> genes;
        std::unique_ptr<GenomeLoader> loader;

        explicit Reference(bool with_neighbours = false) {
            Random random(11);
            std::string fna, map;
            for (int taxid = 1; taxid <= kTaxa; taxid++) {
                for (int gene = 1; gene <= kGenes; gene++) {
                    std::string seq;
                    for (int i = 0; i < kGeneLength; i++) seq += "ACGT"[random.Next(4)];
                    genes[{ taxid, gene }] = seq;
                    std::string const header = ">" + std::to_string(taxid) + "_" + std::to_string(gene) + "\n";
                    size_t const start = fna.size() + header.size();
                    fna += header + seq + "\n";
                    map += std::to_string(taxid) + "\t" + std::to_string(gene) + "\t" + std::to_string(start) + "\t" +
                           std::to_string(start + seq.size()) + "\n";
                }
            }
            auto const fna_path = dir.Write("reference.fna", fna);
            auto const map_path = dir.Write("reference.map", map);
            loader = std::make_unique<GenomeLoader>(fna_path, map_path);
            loader->LoadAllGenomes();
            if (with_neighbours) loader->SetGeneNeighbours(Neighbours());
        }

        std::string Header() const {
            std::string header = "@HD\tVN:1.6\n";
            for (auto const& [key, seq] : genes) {
                header += "@SQ\tSN:" + std::to_string(key.first) + "_" + std::to_string(key.second) + "\tLN:" + std::to_string(seq.size()) + "\n";
            }
            return header;
        }
    };

    std::string Line(std::string const& qname, int flag, std::string const& rname, int pos, int mapq, std::string const& cigar,
                     std::string const& seq, std::string const& tags) {
        return qname + '\t' + std::to_string(flag) + '\t' + rname + '\t' + std::to_string(pos) + '\t' + std::to_string(mapq) + '\t' +
               cigar + "\t*\t0\t0\t" + seq + '\t' + std::string(seq.size(), seq.empty() ? '*' : 'F') + tags + '\n';
    }

    // A record of `length` read bases on gene (taxid, gene) from 0-based `start`: exact, with a mismatch (X),
    // an insertion or a deletion; `misfit`: a mismatch written as M, which does not fit the gene.
    std::string Aligned(Reference const& ref, std::string const& qname, int flag, int taxid, int gene, int start, int length,
                        int kind, int mapq, std::string const& tags, bool misfit = false) {
        auto const& g = ref.genes.at({ taxid, gene });
        std::string seq = g.substr(start, length);
        std::string cigar = std::to_string(length) + "M";
        int const at = length / 2;
        auto other = [](char c) { return c == 'A' ? 'C' : 'A'; };
        if (misfit) {
            seq[at] = other(seq[at]);
        } else if (kind == 1) {
            seq[at] = other(seq[at]);
            cigar = std::to_string(at) + "M1X" + std::to_string(length - at - 1) + "M";
        } else if (kind == 2) {
            seq = g.substr(start, at) + "TT" + g.substr(start + at, length - at - 2);
            cigar = std::to_string(at) + "M2I" + std::to_string(length - at - 2) + "M";
        } else if (kind == 3) {
            seq = g.substr(start, at) + g.substr(start + at + 3, length - at);
            cigar = std::to_string(at) + "M3D" + std::to_string(length - at) + "M";
        }
        return Line(qname, flag, std::to_string(taxid) + "_" + std::to_string(gene), start + 1, mapq, cigar, seq, tags);
    }

    // A SAM of `reads` reads of every kind ProfileSam sees: pairs, single mates, reads with candidates on other
    // taxa (secondary records), long reads in parts (supplementary records), low MAPQ, alternatives (ZA),
    // records of genes the database lacks, records that do not fit their gene (after their taxon has reads),
    // unusable records, blank lines and CRLF line ends.
    std::string Sam(Reference const& ref, int reads, uint64_t seed, bool misfits = true) {
        Random random(seed);
        std::string sam = ref.Header();
        for (int r = 0; r < reads; r++) {
            std::string const name = "read" + std::to_string(r);
            int const taxid = 1 + static_cast<int>(random.Next(kTaxa));
            int const gene = 1 + static_cast<int>(random.Next(kGenes));
            int const start = static_cast<int>(random.Next(kGeneLength - 140));
            int const kind = static_cast<int>(random.Next(4));
            int const mapq = random.Next(10) == 0 ? 2 : 30 + static_cast<int>(random.Next(31));
            std::string const tags = "\tZU:i:" + std::to_string(random.Next(4)) + "\tZT:i:" + std::to_string(random.Next(2)) +
                                     "\tZA:Z:" + (random.Next(3) == 0 ? std::string("*") : std::to_string(1 + random.Next(kTaxa)) + ":" + std::to_string(random.Next(3)));
            switch (random.Next(8)) {
                case 0:  // single-end read
                    sam += Aligned(ref, name, random.Next(2) ? kReverse : 0, taxid, gene, start, 70, kind, mapq, tags);
                    break;
                case 1:  // read 2 alone
                    sam += Aligned(ref, name, kPaired | kRead2, taxid, gene, start, 70, kind, mapq, tags);
                    break;
                case 2: {  // candidates on another taxon, the primary second
                    int const other = 1 + (taxid % kTaxa);
                    sam += Aligned(ref, name, kSecondary, other, gene, start, 70, 0, 0, "\tZU:i:0\tZT:i:0");
                    sam += Aligned(ref, name, 0, taxid, gene, start, 70, kind, mapq, tags);
                    break;
                }
                case 3: {  // a long read in parts on the taxon's genes
                    int const parts = 2 + static_cast<int>(random.Next(3));
                    for (int p = 0; p < parts; p++) {
                        sam += Aligned(ref, name, p == 0 ? 0 : kSupplementary, taxid, 1 + (gene + p) % kGenes, start, 90, kind, mapq, tags);
                    }
                    break;
                }
                case 4:  // a record of a gene the database lacks, or past a gene's end; with a blank line
                    sam += random.Next(2) ? Aligned(ref, name, 0, taxid, gene, start, 70, 0, mapq, tags).replace(name.size() + 3, 1, "9")
                                          : Line(name, 0, "1_1", kGeneLength - 30, 60, "70M", ref.genes.at({ 1, 1 }).substr(0, 70), tags);
                    sam += "\n";
                    break;
                case 5:  // unusable: unmapped, or not a protal gene; CRLF
                    sam += Line(name, 4, "*", 0, 0, "*", "ACGT", "");
                    sam += Line(name + "x", 0, "chr1", 1, 60, "4M", "ACGT", "\r");
                    break;
                default: {  // a pair
                    int const start2 = std::min(start + 40, kGeneLength - 75);
                    bool const misfit = misfits && r > 50 && random.Next(40) == 0;
                    sam += Aligned(ref, name, kPaired | kBothAlign | kRead1 | kMateReverse, taxid, gene, start, 70, kind, mapq, tags, misfit);
                    sam += Aligned(ref, name, kPaired | kBothAlign | kRead2 | kReverse, taxid, gene, start2, 70, (kind + 1) % 4, mapq, tags);
                }
            }
        }
        return sam;
    }

    // Everything a profile holds, in the order of its maps: taxa, their features and genes, the genes' reads,
    // coverage and variants.
    std::string Dump(profiler::MicrobialProfile const& profile) {
        std::ostringstream os;
        for (auto const& [taxid, taxon] : profile.GetTaxa()) {
            os << "taxon " << taxid << " hits " << taxon.TotalHits() << " fragments " << taxon.Fragments() << '\n';
            for (auto const& [name, value] : profiler::TaxonFeatures(taxon)) os << ' ' << name << '=' << profiler::FeatureString(value);
            os << '\n';
            for (auto const& [geneid, gene] : taxon.GetGenes()) {
                os << " gene " << geneid << ' ' << gene.m_mapped_reads << ' ' << gene.m_mapped_length << ' '
                   << profiler::FeatureString(gene.m_identity_bases) << ' ' << profiler::FeatureString(gene.m_ani_sum) << ' '
                   << gene.m_fragments << " identities";
                for (auto const& [identity, length] : gene.m_read_identities) os << ' ' << identity << ':' << length;
                os << "\n  coverage";
                for (auto c : gene.GetStrainLevel().GetSequenceRangeHandler().CalculateCoverageVector2()) os << ' ' << c;
                os << "\n  variants";
                for (auto const& [pos, bin] : gene.GetStrainLevel().GetVariantHandler().GetVariants()) {
                    os << ' ' << pos << ':' << VariantHandler::VariantBinToString(bin);
                }
                os << '\n';
            }
        }
        return os.str();
    }

    struct Profiled {
        std::string error, rejected, log, dump;
        size_t reads = 0, rejected_reads = 0;
    };

    Profiled Profile(Reference const& ref, std::string const& sam, size_t threads, size_t chunk_bytes = 0,
                     std::optional<size_t> decompress_threads = std::nullopt) {
        profiler::Profiler profiler(*ref.loader);
        profiler.SetDepthIdentityMargin(0.08);
        if (chunk_bytes > 0) profiler.SetChunkBytes(chunk_bytes);
        if (decompress_threads) profiler.SetDecompressThreads(*decompress_threads);
        profiler::MicrobialProfile profile(*ref.loader);
        auto genera = std::make_shared<std::vector<uint32_t>>(std::vector<uint32_t>{ 0, 10, 10, 11, 11, 12 });
        profile.SetGenera(genera);
        std::ostringstream rejected;
        testing::internal::CaptureStderr();
        Profiled out;
        out.error = profiler.ProfileSam(sam, profile, std::optional<std::reference_wrapper<std::ostream>>{ rejected }, 2, 2, 0.0, 15, 0,
                                        false, threads);
        out.log = testing::internal::GetCapturedStderr();
        out.rejected = rejected.str();
        out.reads = profiler.Reads();
        out.rejected_reads = profiler.RejectedReads();
        out.dump = Dump(profile);
        return out;
    }

    // The lines of a log that report on the SAM (skipped records, no records, no tags). Records that do not fit
    // their gene are shown only for the first few of the process, so those lines differ between runs.
    std::string Report(std::string const& log) {
        std::istringstream is(log);
        std::string line, report;
        while (std::getline(is, line)) {
            if (line.find("skipped ") != std::string::npos || line.find("contains no usable") != std::string::npos ||
                line.rfind("Warning: ", 0) == 0) {
                report += line + '\n';
            }
        }
        return report;
    }

    void ExpectSame(Profiled const& a, Profiled const& b) {
        EXPECT_EQ(a.error, b.error);
        EXPECT_EQ(a.reads, b.reads);
        EXPECT_EQ(a.rejected_reads, b.rejected_reads);
        EXPECT_EQ(a.rejected, b.rejected);
        EXPECT_EQ(Report(a.log), Report(b.log));
        // After an error the profile is incomplete, and not used.
        if (a.error.empty()) EXPECT_EQ(a.dump, b.dump);
    }

    // The QNAMEs of a text's record lines.
    std::vector<std::string> Names(std::string const& text) {
        std::vector<std::string> names;
        for (size_t at = 0; at < text.size();) {
            size_t end = text.find('\n', at);
            if (end == std::string::npos) end = text.size();
            if (auto name = sam_chunks::RecordName(text.data() + at, text.data() + end)) names.emplace_back(*name);
            at = end + 1;
        }
        return names;
    }
}

TEST(SamChunks, AChunkEndsBeforeARecordOfAnotherRead) {
    using sam_chunks::ChunkEnd;
    std::string const text = "@HD\tVN:1.6\na\t1\nb\t1\nb\t2\n\nc\t1\nc";
    // The last complete record line of another read than the one before it is c's: the chunk ends before it.
    EXPECT_EQ(ChunkEnd(text.data(), text.size()), text.find("c\t1"));
    // Within b's lines, before b.
    std::string const two = "a\t1\nb\t1\nb\t2\n";
    EXPECT_EQ(ChunkEnd(two.data(), two.size()), two.find("b\t1"));
    // One read, or no complete line: no end yet.
    std::string const one = "@HD\nb\t1\nb\t2\r\nb";
    EXPECT_EQ(ChunkEnd(one.data(), one.size()), 0u);
    EXPECT_EQ(ChunkEnd("a\t1", 3), 0u);
    // A line without a tab is named by all of it, CR aside.
    std::string const bare = "a\r\nb\n";
    EXPECT_EQ(ChunkEnd(bare.data(), bare.size()), bare.find('b'));
}

TEST(SamChunks, TheChunksAreTheStreamCutBetweenReads) {
    Reference ref;
    std::string const sam = Sam(ref, 300, 3) + "lastread\t0";  // ends without a newline
    for (size_t bytes : { 1, 7, 100, 1000, 100000 }) {
        SCOPED_TRACE(bytes);
        std::istringstream is(sam);
        sam_chunks::ChunkReader reader(is, bytes, 3);
        std::string joined, last_name;
        size_t lines = 0, index = 0, chunks = 0;
        sam_chunks::Chunk chunk;
        bool last = false;
        while (reader.Next(chunk)) {
            EXPECT_FALSE(last);
            EXPECT_EQ(chunk.index, index++);
            EXPECT_EQ(chunk.first_line, lines);
            // A chunk starts with a read of its own.
            auto const names = Names(chunk.text);
            if (!names.empty()) {
                if (!last_name.empty()) EXPECT_NE(names.front(), last_name);
                last_name = names.back();
            }
            lines += static_cast<size_t>(std::count(chunk.text.begin(), chunk.text.end(), '\n'));
            joined += chunk.text;
            last = chunk.last;
            chunks++;
        }
        if (bytes == 1) EXPECT_GT(chunks, 300u);
        reader.Stop();
        EXPECT_TRUE(last);
        EXPECT_FALSE(reader.Bad());
        EXPECT_EQ(joined, sam);
    }
}

TEST(SamChunks, ParallelForRunsEachOnceAndPassesOnAnException) {
    std::vector<int> runs(1000, 0);
    sam_chunks::ParallelFor(runs.size(), 8, [&](size_t i) { runs[i]++; });
    EXPECT_EQ(std::count(runs.begin(), runs.end(), 1), 1000);
    EXPECT_THROW(sam_chunks::ParallelFor(100, 4, [](size_t i) { if (i == 50) throw std::runtime_error("x"); }), std::runtime_error);
    sam_chunks::ParallelFor(0, 4, [](size_t) { FAIL(); });
}

TEST(ProfileSam, OnSeveralThreadsTheProfileIsTheSame) {
    Reference ref;
    auto const sam = ref.dir.Write("sample.sam", Sam(ref, 3000, 5));
    auto const serial = Profile(ref, sam, 1);
    ASSERT_EQ(serial.error, "");
    EXPECT_GT(serial.rejected_reads, 0u);
    EXPECT_NE(serial.log.find("skipped"), std::string::npos);
    EXPECT_EQ(serial.dump.find("taxon"), 0u);
    for (size_t threads : { 2, 3, 8 }) {
        for (size_t bytes : { 1, 500, 20000, 0 }) {
            SCOPED_TRACE("threads " + std::to_string(threads) + ", chunks of " + std::to_string(bytes) + " bytes");
            ExpectSame(serial, Profile(ref, sam, threads, bytes));
        }
    }
}

TEST(ProfileSam, OnSeveralThreadsTheGeneNeighboursFeaturesAreTheSame) {
    // adjacent_support sums the smoothed shares of the links of a taxon's reads. On several threads they are summed
    // chunk by chunk, on one thread read by read: the sum has to come out the same to the last bit, as the GTDB
    // build's parity check (check_model_parity.py) profiles a training sample alone on all threads, where the
    // collector had profiled it among other samples on fewer.
    Reference ref(true);
    auto const sam = ref.dir.Write("sample.sam", Sam(ref, 3000, 17));
    auto const serial = Profile(ref, sam, 1);
    ASSERT_EQ(serial.error, "");
    size_t judged = 0;
    for (size_t at = serial.dump.find(" adjacent_support="); at != std::string::npos; at = serial.dump.find(" adjacent_support=", at + 1)) {
        judged += serial.dump.compare(at, 22, " adjacent_support=0.5 ") != 0;
    }
    EXPECT_EQ(judged, static_cast<size_t>(kTaxa));  // every taxon has links the clades judge
    for (size_t threads : { 2, 3, 8 }) {
        for (size_t bytes : { 1, 500, 20000, 0 }) {
            SCOPED_TRACE("threads " + std::to_string(threads) + ", chunks of " + std::to_string(bytes) + " bytes");
            ExpectSame(serial, Profile(ref, sam, threads, bytes));
        }
    }
}

TEST(ProfileSam, OnSeveralThreadsTheCompressedSamsGiveTheSameProfile) {
    Reference ref;
    auto const text = Sam(ref, 2000, 9);
    auto const header_end = text.find("\nread0\t") + 1;
    for (auto const* name : { "sample.sam.zst", "sample.sam.gz" }) {
        SCOPED_TRACE(name);
        auto const path = (ref.dir.path / name).string();
        {
            SamOutput out(path, SamCompressionOf(path));
            std::vector<uint64_t> genes;
            out.Write(text.data() + header_end, text.size() - header_end, genes);
            ASSERT_TRUE(out.Finish(text.substr(0, header_end))) << out.Error();
        }
        auto const serial = Profile(ref, path, 1);
        ASSERT_EQ(serial.error, "");
        ExpectSame(serial, Profile(ref, path, 4, 3000));
        ExpectSame(serial, Profile(ref, path, 4));

        // Cut short: the same error.
        fs::resize_file(path, fs::file_size(path) / 2);
        auto const cut = Profile(ref, path, 1);
        EXPECT_NE(cut.error.find("the file is truncated or corrupt"), std::string::npos) << cut.error;
        ExpectSame(cut, Profile(ref, path, 4, 3000));
    }

    // A gzip file of one member (not BGZF), cut in the middle: the stream ends early, which only reading shows.
    auto const path = (ref.dir.path / "plain.sam.gz").string();
    {
        ogzstream os(path.c_str());
        os << text;
    }
    ExpectSame(Profile(ref, path, 1), Profile(ref, path, 4, 3000));
    fs::resize_file(path, fs::file_size(path) / 2);
    auto const cut = Profile(ref, path, 1);
    EXPECT_NE(cut.error.find("the file is truncated or corrupt"), std::string::npos) << cut.error;
    for (size_t bytes : { 1, 3000, 0 }) ExpectSame(cut, Profile(ref, path, 4, bytes));
}

TEST(ProfileSam, AZstdSamOfManyFramesIsDecompressedOnThreadsOfItsOwn) {
    Reference ref;
    auto const text = Sam(ref, 3000, 13);
    auto const header_end = text.find("\nread0\t") + 1;
    auto const path = (ref.dir.path / "frames.sam.zst").string();
    {
        // Many frames: the records written in pieces of whole lines, as the alignment threads hand them over.
        SamOutput out(path, SamCompression::Zstd);
        std::vector<uint64_t> genes;
        for (size_t at = header_end; at < text.size();) {
            size_t end = text.find('\n', std::min(text.size() - 1, at + 9000)) + 1;
            out.Write(text.data() + at, end - at, genes);
            at = end;
        }
        ASSERT_TRUE(out.Finish(text.substr(0, header_end))) << out.Error();
    }
    std::string error;
    ASSERT_GT(zstd::ReadSeekTable(path, error)->frames.size(), 50u);
    auto const serial = Profile(ref, path, 1);
    ASSERT_EQ(serial.error, "");
    for (size_t decompress : { 0, 1, 3 }) {
        SCOPED_TRACE(decompress);
        ExpectSame(serial, Profile(ref, path, 4, 2000, decompress));
        ExpectSame(serial, Profile(ref, path, 2, 0, decompress));
    }

    // A frame in the middle that does not decompress: an error on any number of threads.
    std::string bytes;
    {
        std::ifstream is(path, std::ios::binary);
        bytes.assign(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>());
    }
    auto const table = zstd::ReadSeekTable(path, error);
    auto const& middle = table->frames[table->frames.size() / 2];
    bytes[middle.compressed_offset + middle.compressed_size / 2] ^= 0x5a;
    std::ofstream(path, std::ios::binary) << bytes;
    for (size_t decompress : { 0, 1, 3 }) {
        auto const corrupt = Profile(ref, path, 4, 2000, decompress);
        EXPECT_NE(corrupt.error.find("the file is truncated or corrupt"), std::string::npos) << decompress << ": " << corrupt.error;
    }
    EXPECT_NE(Profile(ref, path, 1).error.find("the file is truncated or corrupt"), std::string::npos);
}

TEST(ProfileSam, ATaxonLeftWithoutRecordsIsProfiledAsOnOneThread) {
    Reference ref;
    // Taxon 4's first record does not fit its gene: on one thread the taxon is removed and made anew by its next
    // record; on several, ProfileSam notices and profiles the file on one thread.
    std::string const tags = "\tZU:i:1\tZT:i:0\tZA:Z:*";
    std::string text = ref.Header();
    text += Aligned(ref, "first", 0, 2, 1, 10, 70, 0, 60, tags);
    text += Aligned(ref, "misfit", 0, 4, 2, 20, 70, 0, 60, tags, true);
    for (int r = 0; r < 200; r++) text += Aligned(ref, "r" + std::to_string(r), 0, 1 + r % kTaxa, 1 + r % kGenes, r % 100, 70, r % 4, 60, tags);
    auto const sam = ref.dir.Write("misfit.sam", text);
    auto const serial = Profile(ref, sam, 1);
    ASSERT_EQ(serial.error, "");
    EXPECT_EQ(serial.rejected_reads, 1u);
    ExpectSame(serial, Profile(ref, sam, 4, 200));
}

TEST(ProfileSam, OnSeveralThreadsTheErrorsAreTheSame) {
    Reference ref;
    auto const good = Sam(ref, 500, 7, false);
    std::string broken = good;
    size_t const at = broken.find("\nread300\t") + 1;
    broken.insert(at, "bad\tline\n");
    std::string other_db = good;
    other_db.insert(other_db.find("@SQ"), "@SQ\tSN:1_1\tLN:999\n");
    for (auto const& [name, text] : std::vector<std::pair<std::string, std::string>>{ { "broken.sam", broken }, { "other.sam", other_db },
                                                                                      { "empty.sam", "" }, { "header.sam", ref.Header() } }) {
        SCOPED_TRACE(name);
        auto const path = ref.dir.Write(name, text);
        auto const serial = Profile(ref, path, 1);
        for (size_t bytes : { 1, 300, 0 }) ExpectSame(serial, Profile(ref, path, 3, bytes));
        if (name == "broken.sam") EXPECT_NE(serial.error.find("line "), std::string::npos) << serial.error;
    }
}
