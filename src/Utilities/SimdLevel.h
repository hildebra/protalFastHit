// SimdLevel.h - which vector instructions protal's own kernels use, chosen at run time.
//
// Three kernels have versions for more than one instruction set: the flex-cell scan of every seed lookup
// (Hash/FlexScan.h), the closed-syncmer scan (SequenceUtils/KmerIterator.h) and 2-bit packing
// (SequenceUtils/PackedSequence.h, AVX2 only). Each uses the highest level the CPU supports, capped by
// PROTAL_SIMD in the environment:
//
//   PROTAL_SIMD=auto     the CPU's highest level (the default; also for an unset or unknown value)
//   PROTAL_SIMD=avx512   at most AVX-512
//   PROTAL_SIMD=avx2     at most AVX2
//   PROTAL_SIMD=scalar   no vector kernels
//
// Every level gives the same outputs; the cap is for comparing runs and for working round a CPU or a
// virtual machine that misreports what it has. The older PROTAL_FLEX_SCAN=scalar still switches off the
// flex scan's vector versions alone.
//
// The AVX-512 level is the set Intel Ice Lake and AMD Zen 4 brought, and later CPUs of both have: AVX-512 F,
// BW, VL, DQ, VBMI, VBMI2 and VPOPCNTDQ. Earlier AVX-512 CPUs (Skylake-X, Cascade Lake) lack the last three
// and use AVX2. libgcc's CPU detection only reports AVX-512 where the operating system saves its registers.
#pragma once

#include <cstdlib>
#include <string_view>

namespace protal::simd {
    enum class Level : int { scalar = 0, avx2 = 1, avx512 = 2 };

    // The target attribute of the AVX-512 kernels (the instruction sets of the AVX-512 level, and those of
    // x86-64-v3 below them).
#define PROTAL_TARGET_AVX512 \
    __attribute__((target("avx2,bmi,bmi2,lzcnt,popcnt,fma,avx512f,avx512bw,avx512vl,avx512dq,avx512vbmi,avx512vbmi2,avx512vpopcntdq")))

#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__))
#define PROTAL_SIMD_X86 1
#endif

    inline char const* Name(Level level) {
        switch (level) {
            case Level::avx512: return "AVX-512";
            case Level::avx2: return "AVX2";
            default: return "scalar";
        }
    }

    // The highest level this CPU (and its operating system) supports.
    inline Level CpuLevel() {
#ifdef PROTAL_SIMD_X86
        __builtin_cpu_init();
        if (!__builtin_cpu_supports("avx2")) return Level::scalar;
        bool const avx512 = __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512bw") &&
                            __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("avx512dq") &&
                            __builtin_cpu_supports("avx512vbmi") && __builtin_cpu_supports("avx512vbmi2") &&
                            __builtin_cpu_supports("avx512vpopcntdq") && __builtin_cpu_supports("bmi2");
        return avx512 ? Level::avx512 : Level::avx2;
#else
        return Level::scalar;
#endif
    }

    // The cap PROTAL_SIMD sets (avx512 when unset, "auto" or unknown).
    inline Level Cap() {
        char const* value = std::getenv("PROTAL_SIMD");
        if (value == nullptr) return Level::avx512;
        std::string_view const v(value);
        if (v == "scalar") return Level::scalar;
        if (v == "avx2") return Level::avx2;
        return Level::avx512;
    }

    inline Level Min(Level a, Level b) { return static_cast<int>(a) < static_cast<int>(b) ? a : b; }

    // The level the kernels use by default: the CPU's, capped by PROTAL_SIMD. Read once.
    inline Level Default() {
        static Level const level = Min(CpuLevel(), Cap());
        return level;
    }

    // `wanted`, lowered to what this CPU supports (for tests and benchmarks that compare levels).
    inline Level Supported(Level wanted) { return Min(wanted, CpuLevel()); }
}
