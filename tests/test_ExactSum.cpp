// ExactSum (Utilities/ExactSum.h): a sum of doubles that is the same in any order of the values, where a double sum
// differs in its last bits.
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>
#include "ExactSum.h"

using protal::ExactSum;

namespace {
    // ANI-like values (0.9 to 1, as CigarANI gives them) and base-weighted identities (up to 15 kb x identity).
    std::vector<double> Values(size_t n, uint32_t seed) {
        std::mt19937 random(seed);
        std::uniform_real_distribution<double> ani(0.9, 1.0), bases(1, 15000);
        std::vector<double> values;
        for (size_t i = 0; i < n; i++) values.push_back(i % 2 ? ani(random) : ani(random) * bases(random));
        return values;
    }

    double DoubleSum(std::vector<double> const& values) {
        double sum = 0;
        for (double v : values) sum += v;
        return sum;
    }

    ExactSum Exact(std::vector<double> const& values) {
        ExactSum sum;
        for (double v : values) sum += v;
        return sum;
    }
}

TEST(ExactSum, TheSumIsTheSameInAnyOrderWhereADoubleSumIsNot) {
    auto values = Values(20000, 3);
    double const exact = Exact(values).Value();
    double const plain = DoubleSum(values);
    bool double_sum_differed = false;
    std::mt19937 random(5);
    for (int round = 0; round < 10; round++) {
        std::shuffle(values.begin(), values.end(), random);
        EXPECT_EQ(Exact(values).Value(), exact);
        double_sum_differed |= DoubleSum(values) != plain;
    }
    std::reverse(values.begin(), values.end());
    EXPECT_EQ(Exact(values).Value(), exact);
    EXPECT_TRUE(double_sum_differed);  // what ExactSum is for
    EXPECT_NEAR(exact, plain, 1e-9 * plain);
}

TEST(ExactSum, PartsAddUpToTheWhole) {
    auto const values = Values(3000, 7);
    ExactSum parts[3];
    for (size_t i = 0; i < values.size(); i++) parts[(i * 7) % 3] += values[i];
    ExactSum whole;
    whole += parts[2];
    whole += parts[0];
    whole += parts[1];
    EXPECT_TRUE(whole == Exact(values));
}

TEST(ExactSum, SmallSumsAreExact) {
    ExactSum sum;
    EXPECT_EQ(sum.Value(), 0.0);
    double plain = 0;
    for (int i = 0; i < 10; i++) {
        sum += 0.1;
        plain += 0.1;
    }
    // Ten times the double 0.1 is 1.0000000000000000555, which rounds to 1; added up in doubles it is 0.9999999999999999.
    EXPECT_EQ(sum.Value(), 1.0);
    EXPECT_NE(plain, 1.0);
    // Less 1, what is left is that excess exactly: 10 x 0.1 is 1 + 2^-54.
    sum += -1.0;
    sum += 2.5;
    sum += -2.5;
    EXPECT_EQ(sum.Value(), 0x1p-54);
    ExactSum large;
    large += 1e15;
    large += 0.5;
    large += 1;
    EXPECT_EQ(large.Value(), 1e15 + 1.5);
}

TEST(ExactSum, BitsBelowTwoToTheMinus64AreDropped) {
    // A value of 2^-11 or more keeps all its bits: 2^-11 + 2^-63 less 2^-11 is 2^-63 exactly.
    ExactSum kept;
    kept += std::nextafter(0x1p-11, 1.0);
    kept += -0x1p-11;
    EXPECT_EQ(kept.Value(), 0x1p-63);
    // Smaller ones lose their bits below 2^-64, each on its own (towards zero), whatever else is in the sum.
    ExactSum small;
    small += 0x1p-65;
    small += 0x1p-65;
    EXPECT_EQ(small.Value(), 0.0);
    small += 0x1.8p-64;  // 1.5 x 2^-64
    EXPECT_EQ(small.Value(), 0x1p-64);
    small += -0x1.8p-64;
    EXPECT_EQ(small.Value(), 0.0);
}

TEST(ExactSum, ValuesUpToTheLimitAddUpExactly) {
    // The largest value below 2^62, twice, and 1: a sum just below 2^63, held exactly and rounded once.
    double const largest = std::nextafter(ExactSum::kMaxMagnitude, 0.0);  // 2^62 - 2^9
    ExactSum sum;
    sum += largest;
    sum += largest;
    sum += 1.0;
    EXPECT_EQ(sum.Value(), 2 * largest);  // 2^63 - 2^10 + 1 is 2^63 - 2^10 as a double
    sum += -largest;
    sum += -largest;
    EXPECT_EQ(sum.Value(), 1.0);
}

// Beyond 2^62, or not finite, a value has no defined conversion to the fixed point: a debug build stops at it.
TEST(ExactSum, AValueBeyondTheLimitStopsADebugBuild) {
#ifdef NDEBUG
    GTEST_SKIP() << "the limit is asserted in debug builds only";
#else
    EXPECT_DEATH({ ExactSum s; s += -ExactSum::kMaxMagnitude; }, "Assertion");
    EXPECT_DEATH({ ExactSum s; s += std::nan(""); }, "Assertion");
#endif
}
