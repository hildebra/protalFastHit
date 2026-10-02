#pragma once

// PROTAL_CLONE_V3 compiles a function twice, for plain x86-64 and for x86-64-v3 (AVX2, BMI1/2, FMA, F16C,
// LZCNT, MOVBE and POPCNT), and the loader runs the copy the CPU supports (GCC's target_clones, through an
// ifunc). This is how one protal binary runs on any x86-64 CPU and still uses the newer instructions where
// they were found to pay (docs/claude/2026-10-02-one-binary): popcount in tsl's sparse maps, which plain
// x86-64 has to compute in software, and the loops the compiler vectorises with AVX2.
//
// - A clone is a call of its own: it is never inlined into its callers. What it inlines is compiled for its
//   target, so it pays on functions that do a lot of work per call, or that would not be inlined anyway.
// - No function-local statics in a cloned function.
// - Floating point rounds the same in both copies: protal is compiled with -ffp-contract=off, so the
//   x86-64-v3 copy does not fuse multiplications and additions into FMA instructions.
//
// Needs GCC 12 or later on x86-64 Linux with glibc (ifunc). Elsewhere, or with -DPROTAL_NO_CLONES, the
// macro is empty and every function is compiled once, for the target of the build.

#include <cstddef>  // defines __GLIBC__ on glibc

#if defined(__x86_64__) && defined(__linux__) && defined(__GLIBC__) && defined(__GNUC__) && !defined(__clang__) && \
    __GNUC__ >= 12 && !defined(PROTAL_NO_CLONES)
#define PROTAL_CLONE_V3 __attribute__((target_clones("default", "arch=x86-64-v3")))
#else
#define PROTAL_CLONE_V3
#endif
