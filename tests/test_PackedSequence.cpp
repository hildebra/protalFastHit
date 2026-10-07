// Unit tests for the 2-bit gene store (SequenceUtils/PackedSequence.h and the genes of GenomeLoader):
// packing and unpacking with the table and with AVX2 give the same bytes and bases, every character
// that is not A, C, G or T is stored as the first base it can stand for, genes filled by several
// threads in pieces equal genes packed whole, and a loader holds its genes packed.
#include <gtest/gtest.h>
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <unistd.h>
#include <vector>
#include "SequenceUtils/GenomeLoader.h"
#include "SequenceUtils/PackedSequence.h"
#include "TestUtil.h"

namespace fs = std::filesystem;
using namespace protal;
using protal::test::ScratchDir;

namespace {
    // What a character becomes: the test's own statement of the rule in PackedSequence.h (an IUPAC code is the first
    // base, in the order A C G T, it stands for; any other character is A).
    char Stored(char c) {
        switch (c) {
            case 'A': case 'a': case 'N': case 'n': case 'R': case 'r': case 'W': case 'w': case 'M': case 'm':
            case 'D': case 'd': case 'H': case 'h': case 'V': case 'v':
                return 'A';
            case 'C': case 'c': case 'Y': case 'y': case 'S': case 's': case 'B': case 'b':
                return 'C';
            case 'G': case 'g': case 'K': case 'k':
                return 'G';
            case 'T': case 't':
                return 'T';
        }
        return 'A';
    }

    std::string StoredOf(std::string const& s) {
        std::string out = s;
        for (auto& c : out) c = Stored(c);
        return out;
    }

    std::string Random(std::mt19937_64& rng, size_t n, std::string const& alphabet) {
        std::string s(n, 'A');
        for (auto& c : s) c = alphabet[rng() % alphabet.size()];
        return s;
    }

    // Both code paths where the CPU has AVX2, else the table twice.
    std::vector<bool> Paths() {
        return packed::CpuHasAvx2() ? std::vector<bool>{ false, true } : std::vector<bool>{ false };
    }

    struct Restore {
        ~Restore() { packed::UseAvx2(true); }
    };
}

TEST(PackedSequence, CodesOfEveryCharacter) {
    for (int c = 0; c < 256; c++) {
        char one = static_cast<char>(c);
        uint8_t byte = 0;
        packed::PackScalar(&one, 1, &byte);
        char back = 0;
        packed::UnpackScalar(&byte, 1, &back);
        EXPECT_EQ(back, Stored(one)) << "character " << c;
    }
}

TEST(PackedSequence, RoundTripOfAllLengthsOnBothPaths) {
    Restore restore;
    std::mt19937_64 rng(1);
    for (bool avx2 : Paths()) {
        packed::UseAvx2(avx2);
        for (size_t n = 0; n <= 700; n++) {
            std::string s = Random(rng, n, "ACGT");
            std::vector<uint8_t> bytes(packed::Bytes(n) + 1, 0xEE);
            packed::Pack(s.data(), n, bytes.data());
            EXPECT_EQ(bytes.back(), 0xEE) << "wrote past the packed bytes, length " << n;
            std::string back(n + 8, '#');
            packed::Unpack(bytes.data(), n, back.data());
            EXPECT_EQ(back.substr(0, n), s) << "length " << n << (avx2 ? " avx2" : " table");
            EXPECT_EQ(back.substr(n), std::string(8, '#')) << "wrote past the bases, length " << n;
        }
        for (size_t n : { 4095u, 4096u, 4097u, 10000u, 65537u }) {
            std::string s = Random(rng, n, "ACGT");
            std::vector<uint8_t> bytes(packed::Bytes(n));
            packed::Pack(s.data(), n, bytes.data());
            std::string back(n, '#');
            packed::Unpack(bytes.data(), n, back.data());
            EXPECT_EQ(back, s) << "length " << n;
        }
    }
}

TEST(PackedSequence, OtherCharactersAreStoredAsTheirFirstBase) {
    Restore restore;
    std::mt19937_64 rng(2);
    for (bool avx2 : Paths()) {
        packed::UseAvx2(avx2);
        // A few ambiguity codes among the bases (the vector path meets them in some blocks), and
        // strings full of them, in both cases.
        for (std::string const alphabet : { std::string("ACGT"), std::string("ACGTN"), std::string("ACGTacgtNRYSWKMBDHVnryswkmbdhv-"),
                                            std::string("NRYSWKMBDHVnryswkmbdhv-") }) {
            for (size_t n : { 1u, 5u, 31u, 32u, 33u, 64u, 127u, 128u, 129u, 500u, 4200u }) {
                for (int rep = 0; rep < 20; rep++) {
                    std::string s = Random(rng, n, alphabet);
                    if (alphabet.size() == 5 && rep % 2) {
                        // mostly bases, one N
                        s = Random(rng, n, "ACGT");
                        s[rng() % n] = 'N';
                    }
                    std::vector<uint8_t> bytes(packed::Bytes(n));
                    packed::Pack(s.data(), n, bytes.data());
                    std::string back(n, '#');
                    packed::Unpack(bytes.data(), n, back.data());
                    ASSERT_EQ(back, StoredOf(s)) << "length " << n << " alphabet " << alphabet << (avx2 ? " avx2" : " table");
                }
            }
        }
    }
}

TEST(PackedSequence, TheTwoPathsGiveTheSameBytes) {
    if (!packed::CpuHasAvx2()) GTEST_SKIP() << "no AVX2";
    Restore restore;
    std::mt19937_64 rng(3);
    for (size_t n : { 0u, 1u, 31u, 32u, 100u, 1000u, 5000u }) {
        for (std::string const alphabet : { std::string("ACGT"), std::string("ACGTacgt"), std::string("ACGTN") }) {
            std::string s = Random(rng, n, alphabet);
            std::vector<uint8_t> table(packed::Bytes(n)), vector(packed::Bytes(n));
            packed::UseAvx2(false);
            packed::Pack(s.data(), n, table.data());
            packed::UseAvx2(true);
            packed::Pack(s.data(), n, vector.data());
            EXPECT_EQ(table, vector) << "length " << n;
        }
    }
}

TEST(PackedSequence, ThreadsFillPiecesOfOneGene) {
    Restore restore;
    std::mt19937_64 rng(4);
    for (bool avx2 : Paths()) {
        packed::UseAvx2(avx2);
        for (size_t n : { 1u, 7u, 100u, 1001u, 7777u }) {
            std::string s = Random(rng, n, "ACGTNRY");
            std::vector<uint8_t> whole(packed::Bytes(n), 0);
            packed::Pack(s.data(), n, whole.data());
            for (int rep = 0; rep < 20; rep++) {
                // Random cut points, so the pieces end inside bytes; one thread per piece.
                std::vector<size_t> cuts{ 0, n };
                for (int k = 0; k < 6; k++) cuts.push_back(rng() % (n + 1));
                std::sort(cuts.begin(), cuts.end());
                std::vector<uint8_t> pieces(packed::Bytes(n), 0);
                std::vector<std::thread> threads;
                for (size_t k = 0; k + 1 < cuts.size(); k++) {
                    threads.emplace_back([&, k] { packed::PackInto(pieces.data(), cuts[k], s.data() + cuts[k], cuts[k + 1] - cuts[k]); });
                }
                for (auto& t : threads) t.join();
                ASSERT_EQ(pieces, whole) << "length " << n << (avx2 ? " avx2" : " table");
            }
        }
    }
}

// A gene of up to kInline bases is decoded into the object itself (on the stack), a longer one into a buffer of its own;
// bases held elsewhere are only viewed.
TEST(GeneSequence, KeepsShortGenesOnTheStackAndLongOnesOnTheHeap) {
    static_assert(!std::is_copy_constructible_v<GeneSequence> && !std::is_move_constructible_v<GeneSequence>,
                  "a view of the sequence points into the object");
    std::mt19937_64 rng(5);
    for (size_t n : std::vector<size_t>{ 0, 1, GeneSequence::kInline, GeneSequence::kInline + 1, 20000 }) {
        std::string s = Random(rng, n, "ACGT");
        std::vector<uint8_t> bytes(packed::Bytes(n) + 1);
        packed::Pack(s.data(), n, bytes.data());
        GeneSequence sequence(bytes.data(), n);
        auto const object = reinterpret_cast<uintptr_t>(&sequence), data = reinterpret_cast<uintptr_t>(sequence.data());
        EXPECT_EQ(data >= object && data + n <= object + sizeof(GeneSequence), n <= GeneSequence::kInline) << "length " << n;
        EXPECT_EQ(sequence.size(), n);
        EXPECT_EQ(sequence.length(), n);
        EXPECT_EQ(std::string_view(sequence), s);
        EXPECT_TRUE(sequence == s) << "length " << n;
        if (n > 10) EXPECT_EQ(sequence.substr(3, 7), std::string_view(s).substr(3, 7));
    }
    std::string external = "ACGTN";
    GeneSequence view{ std::string_view(external) };
    EXPECT_EQ(view.data(), external.data());
}

TEST(PackedSequence, AnyRangeUnpacksAsTheWholeGeneDoes) {
    Restore restore;
    std::mt19937_64 rng(8);
    for (size_t n : { 1u, 5u, 130u, 517u }) {
        std::string const s = Random(rng, n, "ACGT");
        std::vector<uint8_t> bytes(packed::Bytes(n) + 1);
        packed::Pack(s.data(), n, bytes.data());
        for (bool avx2 : Paths()) {
            packed::UseAvx2(avx2);
            for (size_t first = 0; first <= n; first++) {
                for (size_t length : { 0u, 1u, 2u, 3u, 4u, 5u, 31u, 127u, 128u, 129u, 300u }) {
                    if (first + length > n) continue;
                    std::string out(length + 2, '#');  // guard bytes: nothing is written outside the range
                    packed::UnpackRange(bytes.data(), first, length, out.data() + 1);
                    EXPECT_EQ(out.substr(1, length), s.substr(first, length)) << "n " << n << " first " << first << " length " << length;
                    EXPECT_EQ(out.front(), '#');
                    EXPECT_EQ(out.back(), '#');
                }
            }
        }
    }
}

TEST(GeneSequence, AWindowDecodesOnlyItsBasesAtTheirGenePositions) {
    std::mt19937_64 rng(9);
    for (size_t n : { 1u, 200u, 4096u, 4097u, 9000u }) {
        std::string const s = Random(rng, n, "ACGT");
        std::vector<uint8_t> bytes(packed::Bytes(n) + 1);
        packed::Pack(s.data(), n, bytes.data());
        for (auto [begin, end] : std::vector<std::pair<size_t, size_t>>{ { 0, n }, { 0, 1 }, { n / 2, n / 2 + 37 }, { n > 3 ? n - 3 : 0, n },
                                                                           { 7, 3 }, { n, n + 10 }, { 5, 1u << 30 } }) {
            GeneSequence window(bytes.data(), n, begin, end);
            size_t const b = std::min(begin, std::min(end, n)), e = std::min(end, n);
            EXPECT_EQ(window.size(), n);  // as long as the gene, indexed by gene position
            EXPECT_EQ(window.Begin(), b);
            EXPECT_EQ(window.End(), e < b ? b : e);
            for (size_t i = b; i < e; i++) ASSERT_EQ(window[i], s[i]) << "n " << n << " window " << begin << "-" << end << " base " << i;
            if (e > b) EXPECT_EQ(window.substr(b, e - b), std::string_view(s).substr(b, e - b));
        }
    }
    // A whole gene has the window of the whole gene.
    std::string const s = "ACGTACGTAC";
    std::vector<uint8_t> bytes(packed::Bytes(s.size()) + 1);
    packed::Pack(s.data(), s.size(), bytes.data());
    GeneSequence whole(bytes.data(), s.size());
    EXPECT_EQ(whole.Begin(), 0u);
    EXPECT_EQ(whole.End(), s.size());
    EXPECT_EQ(std::string_view(whole), s);
}

TEST(GeneSequence, TheBasesOutsideAWindowAreNotLeftToChance) {
#ifdef PROTAL_GENE_ASAN
    std::mt19937_64 rng(10);
    for (size_t n : { 300u, 5000u }) {  // the stack buffer and the heap one
        std::string const s = Random(rng, n, "ACGT");
        std::vector<uint8_t> bytes(packed::Bytes(n) + 1);
        packed::Pack(s.data(), n, bytes.data());
        char const* data = nullptr;
        {
            GeneSequence window(bytes.data(), n, 100, 200);
            data = window.data();
            EXPECT_FALSE(__asan_address_is_poisoned(data + 100));
            EXPECT_FALSE(__asan_address_is_poisoned(data + 199));
            EXPECT_TRUE(__asan_address_is_poisoned(data + 16));      // before the window (ASan poisons whole 8-byte granules)
            EXPECT_TRUE(__asan_address_is_poisoned(data + n - 16));  // after it
        }
        if (n <= GeneSequence::kInline) {  // on the stack; the poison is gone with the object
            GeneSequence whole(bytes.data(), n);
            EXPECT_FALSE(__asan_address_is_poisoned(whole.data() + 8));
            EXPECT_FALSE(__asan_address_is_poisoned(whole.data() + n - 1));
        }
    }
#else
    GTEST_SKIP() << "only meaningful in an AddressSanitizer build (the undecoded bases are poisoned there)";
#endif
}


namespace {
    // Genes with ambiguity codes and lowercase letters, in the layout of a reference.fna and its map.
    struct AmbiguousReference {
        std::vector<std::string> sequences;
        std::string fna, map;
        AmbiguousReference() {
            std::mt19937_64 rng(6);
            std::ostringstream fna_os, map_os;
            for (int gene = 1; gene <= 40; gene++) {
                std::string s = Random(rng, 3 + rng() % 300, gene % 3 == 0 ? "ACGT" : "ACGTNRYSWKMBDHVacgtn");
                sequences.push_back(s);
                fna_os << '>' << 1 + gene % 3 << '_' << gene << '\n';
                size_t start = fna_os.str().size();
                fna_os << s << '\n';
                map_os << 1 + gene % 3 << '\t' << gene << '\t' << start << '\t' << start + s.size() << '\n';
            }
            fna = fna_os.str();
            map = map_os.str();
        }
    };
}

TEST(GenomeLoaderPacked, GenesAreHeldPackedWithAmbiguityCodesAsTheirFirstBase) {
    ScratchDir dir;
    AmbiguousReference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    for (int threads : { 1, 4 }) {
        GenomeLoader preloaded(fna, map);
        preloaded.LoadAllGenomes(threads);
        GenomeLoader on_demand(fna, map);
        for (int gene = 1; gene <= 40; gene++) {
            int const taxid = 1 + gene % 3;
            std::string const expected = StoredOf(ref.sequences[gene - 1]);
            EXPECT_EQ(preloaded.GeneLength(taxid, gene), expected.size());
            EXPECT_EQ(preloaded.GetGenome(taxid).GetGene(gene).Sequence(), expected) << "gene " << gene << ", " << threads << " threads";
            EXPECT_EQ(on_demand.GetGenome(taxid).GetGeneOMP(gene).Sequence(), expected) << "gene " << gene << " on demand";
        }
    }
}

TEST(GenomeLoaderPacked, AWindowOfAGeneIsTheSameBasesAsTheWholeGene) {
    ScratchDir dir;
    AmbiguousReference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    GenomeLoader preloaded(fna, map);
    preloaded.LoadAllGenomes(2);
    GenomeLoader on_demand(fna, map);  // genes loaded when asked for
    std::mt19937_64 rng(12);
    for (int gene = 1; gene <= 40; gene++) {
        int const taxid = 1 + gene % 3;
        std::string const expected = StoredOf(ref.sequences[gene - 1]);
        for (Gene const* g : { &preloaded.GetGenome(taxid).GetGene(gene), &on_demand.GetGenome(taxid).GetGeneOMP(gene) }) {
            for (int round = 0; round < 8; round++) {
                size_t const begin = rng() % (expected.size() + 1), end = begin + rng() % 120;
                auto const window = g->Window(begin, end);
                ASSERT_EQ(window.size(), expected.size());
                size_t const last = std::min(end, expected.size());
                for (size_t i = begin; i < last; i++) ASSERT_EQ(window[i], expected[i]) << "gene " << gene << " base " << i;
            }
        }
    }
    Gene never_loaded;
    EXPECT_TRUE(never_loaded.Window(0, 10).empty());
}

// The end of --build frees the genes' sequences (GenomeLoader::ReleaseGeneSequences): the genes keep their
// lengths, have no sequence, and reading one again stops protal instead of reading the wrong bytes.
TEST(GenomeLoaderPacked, ReleasedGenesKeepTheirLengthsAndCannotBeReadAgain) {
    ScratchDir dir;
    AmbiguousReference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    GenomeLoader preloaded(fna, map);
    preloaded.LoadAllGenomes(2);
    GenomeLoader on_demand(fna, map);
    ASSERT_FALSE(on_demand.GetGenome(2).GetGeneOMP(1).Sequence().empty());  // taxon 2's genes loaded, 1's and 3's not
    preloaded.ReleaseGeneSequences();
    on_demand.ReleaseGeneSequences();
    for (int gene = 1; gene <= 40; gene++) {
        int const taxid = 1 + gene % 3;
        EXPECT_EQ(preloaded.GeneLength(taxid, gene), ref.sequences[gene - 1].size());
        EXPECT_FALSE(preloaded.GetGenome(taxid).GetGene(gene).IsLoaded());
        EXPECT_TRUE(preloaded.GetGenome(taxid).GetGeneOMP(gene).Sequence().empty()) << gene;
    }
    EXPECT_TRUE(on_demand.GetGenome(2).GetGeneOMP(1).Sequence().empty());
    EXPECT_EXIT(on_demand.GetGenome(3).GetGeneOMP(2), testing::ExitedWithCode(EX_SOFTWARE), "freed at the end of the build");
}

// The flat tables (GetGenome, GetGeneOMP(taxid, gene), PreloadedGene) reach the very genomes and genes the genome map
// holds: preloaded (the gene table), loaded on demand (through the genome, which is loaded on first use), and after the
// sequences are freed; a gene id the genome does not have fails as through the genome, an unknown taxid stops protal.
TEST(GenomeLoaderPacked, TheFlatTablesReachTheGenomesAndGenesOfTheMap) {
    ScratchDir dir;
    AmbiguousReference ref;
    auto fna = dir.Write("reference.fna", ref.fna);
    auto map = dir.Write("reference.map", ref.map);
    GenomeLoader preloaded(fna, map);
    GenomeLoader on_demand(fna, map);
    EXPECT_EQ(preloaded.PreloadedGene(2, 1), nullptr) << "no gene table before the preload";
    preloaded.LoadAllGenomes(2);
    for (int gene = 1; gene <= 40; gene++) {
        int const taxid = 1 + gene % 3;
        std::string const expected = StoredOf(ref.sequences[gene - 1]);
        Genome& genome = preloaded.GetGenomeMap().find(taxid).value();
        EXPECT_EQ(&preloaded.GetGenome(taxid), &genome);
        EXPECT_EQ(&static_cast<GenomeLoader const&>(preloaded).GetGenome(taxid), &genome);
        EXPECT_EQ(&preloaded.GetGeneOMP(taxid, gene), &genome.GetGene(gene)) << gene;
        EXPECT_EQ(preloaded.PreloadedGene(taxid, gene), &genome.GetGene(gene)) << gene;
        EXPECT_EQ(preloaded.GetGeneOMP(taxid, gene).Sequence(), expected) << gene;
        EXPECT_EQ(on_demand.PreloadedGene(taxid, gene), nullptr);
        EXPECT_EQ(on_demand.GetGeneOMP(taxid, gene).Sequence(), expected) << gene << " on demand";
        EXPECT_EQ(&on_demand.GetGeneOMP(taxid, gene), &on_demand.GetGenome(taxid).GetGene(gene));
        // Prefetching is harmless for any gene, known or not.
        preloaded.PrefetchGene(taxid, gene);
        preloaded.GetGeneOMP(taxid, gene).PrefetchBases(0, 1000);
    }
    preloaded.PrefetchGene(99, 1);
    preloaded.PrefetchGene(2, 0);
    EXPECT_EQ(preloaded.PreloadedGene(2, 0), nullptr);
    EXPECT_EQ(preloaded.PreloadedGene(1u << 21, 1), nullptr);
    // Gene 0 and a gene past the genome's list: what Genome::GetGeneOMP does (its gene list's at()).
    EXPECT_THROW(preloaded.GetGeneOMP(2, 0), std::out_of_range);
    EXPECT_THROW(preloaded.GetGeneOMP(2, 1000), std::out_of_range);
    EXPECT_EXIT(preloaded.GetGeneOMP(99, 1), testing::ExitedWithCode(10), "");
    EXPECT_EXIT(preloaded.GetGenome(1u << 21), testing::ExitedWithCode(10), "");
    preloaded.ReleaseGeneSequences();
    EXPECT_TRUE(preloaded.GetGeneOMP(2, 1).Sequence().empty());
    EXPECT_EQ(preloaded.GetGeneOMP(2, 1).GetLength(), ref.sequences[0].size());
}

TEST(GenomeLoaderPacked, GenesCutByFramesAndFilledByManyThreadsEqualWholeGenes) {
    // A seekable zstd reference in frames of 1000 bytes, which cut most genes: the threads fill the
    // pieces of a gene concurrently. The genes must come out as when one thread reads a raw file,
    // lower case and ambiguity codes stored as their (first) base.
    ScratchDir dir;
    std::mt19937_64 rng(7);
    std::ostringstream fna_os, map_os;
    std::vector<std::string> sequences;
    for (int gene = 1; gene <= 3000; gene++) {
        std::string s = Random(rng, 500 + rng() % 1500, "ACGT");
        if (gene % 10 == 0) s[rng() % s.size()] = 'N';
        if (gene % 7 == 0) s[rng() % s.size()] = 'y';
        if (gene % 3 == 0) s[rng() % s.size()] = "acgt"[rng() % 4];
        sequences.push_back(s);
        fna_os << '>' << 1 + gene % 50 << '_' << 1 + gene / 50 << '\n';
        size_t start = fna_os.str().size();
        fna_os << s << '\n';
        map_os << 1 + gene % 50 << '\t' << 1 + gene / 50 << '\t' << start << '\t' << start + s.size() << '\n';
    }
    auto fna = dir.Write("reference.fna", fna_os.str());
    auto map = dir.Write("reference.map", map_os.str());
    std::string error;
    ASSERT_TRUE(zstd::CompressFile(fna, fna + ".zst", { 3, 0, 4, 1000 }, true, error)) << error;
    GenomeLoader raw(fna, map), cut(fna + ".zst", map);
    raw.LoadAllGenomes(1);
    cut.LoadAllGenomes(8);
    for (int gene = 1; gene <= 3000; gene++) {
        int const taxid = 1 + gene % 50, id = 1 + gene / 50;
        EXPECT_EQ(cut.GetGenome(taxid).GetGene(id).Sequence(), StoredOf(sequences[gene - 1])) << gene;
        EXPECT_EQ(raw.GetGenome(taxid).GetGene(id).Sequence(), StoredOf(sequences[gene - 1])) << gene;
    }
}
