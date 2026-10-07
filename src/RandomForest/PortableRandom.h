// SPDX-License-Identifier: GPL-2.0-only
#pragma once

// Random draws that are the same with every C++ standard library, for the community designs of simulate_metagenomes.
// The numbers of std::mt19937_64 are fixed by the standard, but its distributions (std::shuffle,
// std::uniform_int_distribution, std::lognormal_distribution, std::poisson_distribution, ...) are each library's own,
// so a seed gave other communities with another library or version. These take the engine's 64-bit numbers alone
// (any engine with 64-bit results: std::mt19937_64, LongRng). Only libm's exp, log and lgamma remain, whose last bits can
// differ between platforms.

#include <cmath>
#include <cstdint>
#include <numbers>
#include <utility>
#include <vector>

namespace protal::sim::portable {

// [0, 1), 53 bits.
template <typename Rng>
double Uniform(Rng& rng) {
    return static_cast<double>(static_cast<std::uint64_t>(rng()) >> 11) * 0x1.0p-53;
}

// [0, n), n > 0, unbiased (Lemire's multiply and reject).
template <typename Rng>
std::uint64_t Below(Rng& rng, std::uint64_t n) {
    __uint128_t m = static_cast<__uint128_t>(static_cast<std::uint64_t>(rng())) * n;
    auto low = static_cast<std::uint64_t>(m);
    if (low < n) {
        std::uint64_t const threshold = -n % n;
        while (low < threshold) {
            m = static_cast<__uint128_t>(static_cast<std::uint64_t>(rng())) * n;
            low = static_cast<std::uint64_t>(m);
        }
    }
    return static_cast<std::uint64_t>(m >> 64);
}

// [lo, hi], lo <= hi.
template <typename Rng>
std::uint64_t Between(Rng& rng, std::uint64_t lo, std::uint64_t hi) {
    return lo + Below(rng, hi - lo + 1);
}

// Fisher-Yates, from the back.
template <typename Rng, typename T>
void Shuffle(std::vector<T>& values, Rng& rng) {
    for (std::size_t i = values.size(); i > 1; --i) {
        std::swap(values[i - 1], values[Below(rng, i)]);
    }
}

// N(0, 1) by Box-Muller, one number of two uniform ones (no state kept between calls).
template <typename Rng>
double Normal(Rng& rng) {
    double const u = 1.0 - Uniform(rng);  // (0, 1]
    double const v = Uniform(rng);
    return std::sqrt(-2.0 * std::log(u)) * std::cos(2.0 * std::numbers::pi * v);
}

// Gamma(shape, scale) by Marsaglia and Tsang.
template <typename Rng>
double Gamma(Rng& rng, double shape, double scale) {
    if (shape < 1.0) {
        double const u = 1.0 - Uniform(rng);
        return Gamma(rng, shape + 1.0, scale) * std::pow(u, 1.0 / shape);
    }
    double const d = shape - 1.0 / 3.0, c = 1.0 / std::sqrt(9.0 * d);
    for (;;) {
        double x, v;
        do {
            x = Normal(rng);
            v = 1.0 + c * x;
        } while (v <= 0.0);
        v = v * v * v;
        double const u = Uniform(rng);
        if (u < 1.0 - 0.0331 * x * x * x * x) return d * v * scale;
        if (u > 0 && std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v))) return d * v * scale;
    }
}

// Poisson(lambda): multiplication below 10, Hoermann's transformed rejection (PTRS, 1993) above, as numpy draws it.
template <typename Rng>
std::uint64_t Poisson(Rng& rng, double lambda) {
    if (!(lambda > 0)) return 0;
    if (lambda < 10) {
        double const limit = std::exp(-lambda);
        double product = 1.0;
        std::uint64_t k = 0;
        for (;;) {
            product *= Uniform(rng);
            if (product <= limit) return k;
            ++k;
        }
    }
    double const slam = std::sqrt(lambda), loglam = std::log(lambda);
    double const b = 0.931 + 2.53 * slam, a = -0.059 + 0.02483 * b;
    double const invalpha = 1.1239 + 1.1328 / (b - 3.4), vr = 0.9277 - 3.6224 / (b - 2);
    for (;;) {
        double const u = Uniform(rng) - 0.5, v = Uniform(rng);
        double const us = 0.5 - std::fabs(u);
        double const k = std::floor((2 * a / us + b) * u + lambda + 0.43);
        if (us >= 0.07 && v <= vr) return static_cast<std::uint64_t>(k);
        if (k < 0 || (us < 0.013 && v > us)) continue;
        if (std::log(v) + std::log(invalpha) - std::log(a / (us * us) + b) <= -lambda + k * loglam - std::lgamma(k + 1)) {
            return static_cast<std::uint64_t>(k);
        }
    }
}

// Negative binomial: the failures before r successes of probability p, as a Poisson of a gamma mean.
template <typename Rng>
std::uint64_t NegativeBinomial(Rng& rng, double r, double p) {
    return Poisson(rng, Gamma(rng, r, (1.0 - p) / p));
}

}  // namespace protal::sim::portable
