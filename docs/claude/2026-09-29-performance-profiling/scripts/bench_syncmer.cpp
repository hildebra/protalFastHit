// Syncmer extraction: protal's SimpleKmerHandler<ClosedSyncmer> against a version that computes each
// read's 7-mers once per strand. Checks the k-mer lists are identical read by read and times both.
//   bench_syncmer READS.fq [READS2.fq ...]
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include "SequenceUtils/KmerIterator.h"

using namespace protal;

// Same syncmers as SimpleKmerHandler{31, 15, ClosedSyncmer{15, 7, 2, full}}: for each 31-mer window,
// the canonical 15-mer core (the middle 15 bases of the 31-mer or of its reverse complement,
// whichever is smaller; the reverse on a tie) is a closed syncmer if the first minimum of its nine
// 7-mers is the 3rd or the 7th.
class FastSyncmers {
    static constexpr int k = 31, s = 7;
    static constexpr uint64_t kmask = (uint64_t{1} << (2 * k)) - 1, cmask = (uint64_t{1} << 30) - 1;
    std::array<uint8_t, 256> m_fwd{}, m_rev{};
    std::vector<uint16_t> m_f7, m_r7;
    uint32_t m_smask;

public:
    explicit FastSyncmers(bool full_smer_mask) : m_smask(ClosedSyncmer::SmerMask(s, full_smer_mask)) {
        for (int c = 0; c < 256; c++) {
            m_fwd[c] = static_cast<uint8_t>(KmerUtils::BaseToInt(static_cast<char>(c), 0));
            m_rev[c] = static_cast<uint8_t>(KmerUtils::BaseToIntC(static_cast<char>(c), 0));
        }
    }

    void operator()(std::string_view seq, KmerList& list) {
        list.clear();
        size_t const n = seq.size();
        if (n < static_cast<size_t>(k)) return;
        // f7[j]: the 7-mer at j, first base in the high bits; r7[j]: its reverse complement, encoded
        // as the reverse k-mer is (so that N reads as 0 on both strands, as in SimpleKmerHandler).
        m_f7.resize(n - s + 1);
        m_r7.resize(n - s + 1);
        uint32_t f = 0, r = 0;
        for (size_t i = 0; i < n; i++) {
            auto const c = static_cast<unsigned char>(seq[i]);
            f = ((f << 2) | m_fwd[c]) & 0x3FFF;
            r = (r >> 2) | (uint32_t{m_rev[c]} << 12);
            if (i + 1 >= static_cast<size_t>(s)) {
                m_f7[i + 1 - s] = static_cast<uint16_t>(f & m_smask);
                m_r7[i + 1 - s] = static_cast<uint16_t>(r & m_smask);
            }
        }
        uint64_t kf = 0, kr = 0;
        for (size_t i = 0; i < n; i++) {
            auto const c = static_cast<unsigned char>(seq[i]);
            kf = ((kf << 2) | m_fwd[c]) & kmask;
            kr = (kr >> 2) | (uint64_t{m_rev[c]} << (2 * (k - 1)));
            if (i + 1 < static_cast<size_t>(k)) continue;
            size_t const p = i + 1 - k;  // window start
            uint64_t const core_f = (kf >> 16) & cmask, core_r = (kr >> 16) & cmask;
            bool const forward = core_f < core_r;
            // The core's 7-mers in its own order: forward f7[p+8+j]; reverse r7[p+16-j].
            uint16_t const* v = forward ? &m_f7[p + 8] : &m_r7[p + 8];
            int best = 0;
            uint16_t min = forward ? v[0] : v[8];
            for (int j = 1; j < 9; j++) {
                uint16_t const x = forward ? v[j] : v[8 - j];
                if (x < min) { min = x; best = j; }
            }
            if (best == 2 || best == 6) list.emplace_back(forward ? kf : kr, p);
        }
    }
};

// As FastSyncmers, without data-dependent branches: each 7-mer is packed with its index in the core
// (value << 4 | index), so the smallest packed value is the first minimum; both strands are
// evaluated and the canonical one selected.
class BranchFreeSyncmers {
    static constexpr int k = 31, s = 7;
    static constexpr uint64_t kmask = (uint64_t{1} << (2 * k)) - 1, cmask = (uint64_t{1} << 30) - 1;
    std::array<uint8_t, 256> m_fwd{}, m_rev{};
    std::vector<uint32_t> m_f7, m_r7;
    uint32_t m_smask;

public:
    explicit BranchFreeSyncmers(bool full_smer_mask) : m_smask(ClosedSyncmer::SmerMask(s, full_smer_mask)) {
        for (int c = 0; c < 256; c++) {
            m_fwd[c] = static_cast<uint8_t>(KmerUtils::BaseToInt(static_cast<char>(c), 0));
            m_rev[c] = static_cast<uint8_t>(KmerUtils::BaseToIntC(static_cast<char>(c), 0));
        }
    }

    void operator()(std::string_view seq, KmerList& list) {
        list.clear();
        size_t const n = seq.size();
        if (n < static_cast<size_t>(k)) return;
        m_f7.resize(n - s + 1);
        m_r7.resize(n - s + 1);
        uint32_t f = 0, r = 0;
        for (size_t i = 0; i < n; i++) {
            auto const c = static_cast<unsigned char>(seq[i]);
            f = ((f << 2) | m_fwd[c]) & 0x3FFF;
            r = (r >> 2) | (uint32_t{m_rev[c]} << 12);
            if (i + 1 >= static_cast<size_t>(s)) {
                m_f7[i + 1 - s] = (f & m_smask) << 4;
                m_r7[i + 1 - s] = (r & m_smask) << 4;
            }
        }
        uint64_t kf = 0, kr = 0;
        for (size_t i = 0; i < n; i++) {
            auto const c = static_cast<unsigned char>(seq[i]);
            kf = ((kf << 2) | m_fwd[c]) & kmask;
            kr = (kr >> 2) | (uint64_t{m_rev[c]} << (2 * (k - 1)));
            if (i + 1 < static_cast<size_t>(k)) continue;
            size_t const p = i + 1 - k;
            uint64_t const core_f = (kf >> 16) & cmask, core_r = (kr >> 16) & cmask;
            uint32_t const* fv = &m_f7[p + 8];
            uint32_t const* rv = &m_r7[p + 8];
            uint32_t fmin = fv[0], rmin = rv[8];
            for (uint32_t j = 1; j < 9; j++) {
                fmin = std::min(fmin, fv[j] | j);
                rmin = std::min(rmin, rv[8 - j] | j);
            }
            uint32_t const best = (core_f < core_r ? fmin : rmin) & 0xF;
            if (best == 2 || best == 6) list.emplace_back(core_f < core_r ? kf : kr, p);
        }
    }
};

int main(int argc, char** argv) {
    std::vector<std::string> reads;
    for (int a = 1; a < argc; a++) {
        std::ifstream in(argv[a]);
        std::string line;
        for (size_t l = 0; std::getline(in, line); l++) if (l % 4 == 1) reads.push_back(line);
    }
    std::printf("%zu reads\n", reads.size());

    ClosedSyncmer minimizer{15, 7, 2, true};
    SimpleKmerHandler<ClosedSyncmer> current{31, 15, minimizer};
    FastSyncmers fast{true};
    BranchFreeSyncmers branchfree{true};
    KmerList a, b, c;

    // Equality, read by read (and the reads reversed-complemented, and with Ns).
    size_t mismatched = 0, mismatched_bf = 0, total = 0;
    for (auto const& read : reads) {
        for (int variant = 0; variant < 3; variant++) {
            std::string seq = read;
            if (variant == 1) seq = KmerUtils::ReverseComplement(seq);
            if (variant == 2) for (size_t i = 0; i < seq.size(); i += 37) seq[i] = 'N';
            current(std::string_view(seq), a);
            fast(std::string_view(seq), b);
            branchfree(std::string_view(seq), c);
            mismatched += a != b;
            mismatched_bf += a != c;
            total += a.size();
        }
    }
    std::printf("syncmers %zu, reads with different lists: fast %zu, branch-free %zu\n", total, mismatched, mismatched_bf);

    auto time = [&](auto& handler, KmerList& list) {
        size_t sum = 0;
        auto t0 = std::chrono::steady_clock::now();
        for (auto const& read : reads) { handler(std::string_view(read), list); sum += list.size(); }
        return std::pair{ std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), sum };
    };
    for (int rep = 0; rep < 3; rep++) {
        auto [cur, s1] = time(current, a);
        auto [fst, s2] = time(fast, b);
        auto [bf, s3] = time(branchfree, c);
        std::printf("ns/read: current %.0f, fast %.0f (%.1fx), branch-free %.0f (%.1fx)  [%zu %zu %zu]\n",
                    1e9 * cur / reads.size(), 1e9 * fst / reads.size(), cur / fst, 1e9 * bf / reads.size(), cur / bf, s1, s2, s3);
    }
}
