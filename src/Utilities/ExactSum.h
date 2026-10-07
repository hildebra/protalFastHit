#pragma once

// A sum of doubles that comes out the same in whatever order they are added. A floating-point sum is rounded after
// every addition, so its last bits depend on the order of the values: the records of a SAM that multi-threaded alignment
// wrote come in a different order from run to run, and a sum over them (an ANI sum, say) changed in its last printed
// digit between runs of the same build. ExactSum holds each value in fixed point, in units of 2^-64, in a 128-bit
// integer, and integers add up to the same in any order. Value() is the exact sum, rounded once.
//
// Exact for the values protal sums (identities and ANIs, bases times identity): a value of 2^-11 or more has no bits
// below 2^-64; the bits below 2^-64 of a smaller one are dropped, the same in any order. Values must be finite and below
// 2^62 in magnitude, and so must the sum (2^63). Debug builds assert the values' limit: beyond it, the conversion to
// __int128 is undefined.

#include <cassert>
#include <cmath>

namespace protal {
    class ExactSum {
        static constexpr double kUnit = 18446744073709551616.0;  // 2^64
        __int128 m_units = 0;

    public:
        static constexpr double kMaxMagnitude = 4611686018427387904.0;  // 2^62: values are below it

        ExactSum& operator+=(double value) {
            assert(std::isfinite(value) && std::fabs(value) < kMaxMagnitude);
            m_units += static_cast<__int128>(value * kUnit);
            return *this;
        }

        ExactSum& operator+=(ExactSum const& other) {
            m_units += other.m_units;
            return *this;
        }

        double Value() const {
            return static_cast<double>(m_units) / kUnit;
        }

        bool operator==(ExactSum const& other) const {
            return m_units == other.m_units;
        }
    };
}
