// Unit tests for the reverse complement of reads (KmerUtils::ReverseComplement, ...Into): the
// definition (ACGT swap, every other byte becomes N, order reversed) for every byte and for random
// reads, and the buffer the Into form writes into.
#include <gtest/gtest.h>
#include <random>
#include <string>
#include <string_view>
#include "SequenceUtils/KmerUtils.h"

namespace {
    // The definition, one base at a time (what KmerUtils::ReverseComplement did before the table).
    std::string Naive(std::string const& forward) {
        std::string reverse;
        for (int i = static_cast<int>(forward.size()) - 1; i >= 0; i--) reverse += KmerUtils::Complement(forward[i]);
        return reverse;
    }
}

TEST(ReverseComplement, EveryByteHasItsComplement) {
    for (int b = 0; b < 256; b++) {
        std::string const one(1, static_cast<char>(b));
        std::string const expected = (b == 'A' ? "T" : b == 'C' ? "G" : b == 'G' ? "C" : b == 'T' ? "A" : "N");
        EXPECT_EQ(KmerUtils::ReverseComplement(one), expected) << "byte " << b;
    }
}

TEST(ReverseComplement, RandomReadsOfAnyLengthMatchTheDefinition) {
    std::mt19937 rng(3);
    std::string const alphabet = "ACGTACGTACGTNnacgtRYKM-.";
    for (size_t length = 0; length < 400; length++) {
        std::string read(length, 'A');
        for (auto& c : read) c = alphabet[rng() % alphabet.size()];
        EXPECT_EQ(KmerUtils::ReverseComplement(read), Naive(read)) << "length " << length;
        EXPECT_EQ(KmerUtils::ReverseComplement(std::string_view(read)), Naive(read));
    }
}

TEST(ReverseComplement, IntoReusesTheBufferAndResizesIt) {
    std::string buffer;
    KmerUtils::ReverseComplementInto("AACCGGTT", buffer);
    EXPECT_EQ(buffer, "AACCGGTT");
    KmerUtils::ReverseComplementInto(std::string(200, 'A'), buffer);
    EXPECT_EQ(buffer, std::string(200, 'T'));
    char const* const data = buffer.data();
    KmerUtils::ReverseComplementInto("GATTACA", buffer);  // shorter: the same allocation, the right length
    EXPECT_EQ(buffer, "TGTAATC");
    EXPECT_EQ(buffer.data(), data);
    KmerUtils::ReverseComplementInto("", buffer);
    EXPECT_TRUE(buffer.empty());
}
