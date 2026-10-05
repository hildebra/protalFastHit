// Sorting a read's seeds (LookupResult, 16 bytes) by (taxid, gene, read position): std::sort with the comparator (now),
// and by a 128-bit key holding every field (a total order), with std::sort, with insertion sort for small lists, and an
// LSD radix sort. Lists like a paired-end mate's at r226: ~43 seeds on a few taxa and genes, 22 read positions.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>
struct LookupResult { uint32_t taxid, geneid, genepos; uint16_t readpos; bool unique, unique_dist_two; };
static bool Cmp(LookupResult const& a, LookupResult const& b) {
    if (a.taxid != b.taxid) return a.taxid < b.taxid;
    if (a.geneid != b.geneid) return a.geneid < b.geneid;
    return a.readpos < b.readpos;
}
using u128 = unsigned __int128;
static inline u128 Key(LookupResult const& s) {
    return (u128(s.taxid) << 58) | (u128(s.geneid) << 38) | (u128(s.readpos) << 22) | (u128(s.genepos) << 2) | (u128(s.unique) << 1) | u128(s.unique_dist_two);
}
static inline LookupResult Decode(u128 k) {
    LookupResult s; s.taxid = uint32_t(k >> 58); s.geneid = uint32_t(k >> 38) & 0xFFFFF; s.readpos = uint16_t(k >> 22);
    s.genepos = uint32_t(k >> 2) & 0xFFFFF; s.unique = (k >> 1) & 1; s.unique_dist_two = k & 1; return s;
}
static void KeySort(std::vector<LookupResult>& v, std::vector<u128>& keys) {
    keys.resize(v.size()); for (size_t i = 0; i < v.size(); i++) keys[i] = Key(v[i]);
    std::sort(keys.begin(), keys.end()); for (size_t i = 0; i < v.size(); i++) v[i] = Decode(keys[i]);
}
static void KeyInsertion(std::vector<LookupResult>& v, std::vector<u128>& keys) {
    size_t const n = v.size(); keys.resize(n);
    for (size_t i = 0; i < n; i++) { u128 k = Key(v[i]); size_t j = i; while (j > 0 && keys[j - 1] > k) { keys[j] = keys[j - 1]; j--; } keys[j] = k; }
    for (size_t i = 0; i < n; i++) v[i] = Decode(keys[i]);
}
static void KeyHybrid(std::vector<LookupResult>& v, std::vector<u128>& keys) { if (v.size() <= 48) KeyInsertion(v, keys); else KeySort(v, keys); }
// Seeds as a stable sort of the struct by the full order: std::sort of the struct with a total-order comparator.
static bool CmpTotal(LookupResult const& a, LookupResult const& b) { return Key(a) < Key(b); }
int main() {
    std::mt19937 rng(5);
    std::vector<std::vector<LookupResult>> lists(200000);
    for (auto& l : lists) {
        std::geometric_distribution<int> geo(1.0 / 43); size_t n = std::min(256, 1 + geo(rng));
        uint32_t taxa[6]; for (auto& t : taxa) t = rng() % 140000;
        for (size_t i = 0; i < n; i++) {
            LookupResult s; s.taxid = taxa[std::min<uint32_t>(rng() % 8, 5)]; s.geneid = 1 + rng() % 3; s.readpos = uint16_t((rng() % 22) * 6);
            s.genepos = 300 + s.readpos + (rng() % 10 == 0 ? rng() % 900 : 0); s.unique = rng() & 1; s.unique_dist_two = rng() & 1; l.push_back(s);
        }
    }
    size_t total = 0; for (auto& l : lists) total += l.size();
    std::vector<u128> keys;
    auto run = [&](char const* what, auto&& sort) {
        double best = 1e9;
        for (int rep = 0; rep < 5; rep++) { auto copy = lists; auto t = std::chrono::steady_clock::now();
            for (auto& l : copy) sort(l); best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count()); }
        printf("%-44s %.1f ns per seed\n", what, best * 1e9 / total);
    };
    printf("%zu lists, %.1f seeds each\n", lists.size(), double(total) / lists.size());
    run("std::sort, comparator (now)", [&](auto& l) { std::sort(l.begin(), l.end(), Cmp); });
    run("std::sort, total-order comparator", [&](auto& l) { std::sort(l.begin(), l.end(), CmpTotal); });
    run("128-bit keys, std::sort", [&](auto& l) { KeySort(l, keys); });
    run("128-bit keys, insertion sort", [&](auto& l) { KeyInsertion(l, keys); });
    run("128-bit keys, insertion to 48, else std::sort", [&](auto& l) { KeyHybrid(l, keys); });
    // The total order's sorts agree.
    auto a = lists, b = lists, c = lists; for (auto& l : a) KeySort(l, keys); for (auto& l : b) KeyHybrid(l, keys);
    for (auto& l : c) std::sort(l.begin(), l.end(), CmpTotal);
    bool same = true; for (size_t i = 0; i < a.size(); i++) for (size_t j = 0; j < a[i].size(); j++) same &= Key(a[i][j]) == Key(b[i][j]) && Key(a[i][j]) == Key(c[i][j]);
    printf("total-order sorts agree: %s\n", same ? "yes" : "NO");
}
