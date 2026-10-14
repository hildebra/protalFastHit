#!/bin/bash
# Variants of the clone test: the class template's cloned member called from an unused function through a lambda, an
# OpenMP loop, a std::function, or a static member function.
T=~/ctime/clonetest
mkdir -p $T
cd $T
for v in lambda omp function nonclone_omp; do
cat > u_$v.cpp <<EOF
#include <cstddef>
#include <functional>
#include <vector>
#define CLONE __attribute__((target_clones("default", "arch=x86-64-v3")))
inline double Helper(double const* x, size_t n) { double s = 0; for (size_t i = 0; i < n; ++i) s += x[i] * x[i]; return s; }
template<typename T> struct U {
#if defined(V_nonclone_omp)
    T Member(T const* x, size_t n) const { return Helper(x, n) + 4; }
#else
    CLONE T Member(T const* x, size_t n) const { return Helper(x, n) + 4; }
#endif
    T Run(std::vector<T> const& v) const {
#if defined(V_lambda)
        auto f = [&](size_t i) { return Member(v.data() + i, 1); };
        T s = 0; for (size_t i = 0; i < v.size(); ++i) s += f(i); return s;
#elif defined(V_function)
        std::function<T(size_t)> f = [&](size_t i) { return Member(v.data() + i, 1); };
        T s = 0; for (size_t i = 0; i < v.size(); ++i) s += f(i); return s;
#else
        T s = 0;
        #pragma omp parallel for reduction(+:s)
        for (size_t i = 0; i < v.size(); ++i) s += Member(v.data() + i, 1);
        return s;
#endif
    }
};
static double Unused(std::vector<double> const& v) { U<double> u; return u.Run(v); }
int main() { return 0; }
EOF
g++ -O3 -std=c++20 -fopenmp -march=x86-64 -DV_$v -c u_$v.cpp -o u_$v.o
echo "-- $v: $(nm -C u_$v.o | grep -cE 'Helper|Member|Run') symbols: $(nm -C u_$v.o | grep -E 'Helper|Member|Run' | awk '{$1=""; print}' | tr '\n' ';' | cut -c1-300)"
done
