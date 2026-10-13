// Stage1Probe.h - what the first stage of a two-stage protal would keep (--stage1_probe), measured in an ordinary run.
//
// A two-stage protal would look a sample's reads up in a sparse index first: the cores (the 15 bases a k-mer is keyed
// by) whose hash falls in 1/c of its range, each with all its entries. A read would go on to the full alignment, against
// the taxa the first stage found, only if at least m of its lookups there match closely: their best flex cell shares at
// least T of the 16 flank bases (docs/claude/2026-10-12-profiler-fork-ideas.md, follow-up). The full index holds the
// sparse one's entries of every core in it, so its lookups give the first stage's scores: with --stage1_probe, every
// lookup's best score and its core's level (Level: in the sparse index of 1/c for c = 2^level and every smaller c) are
// counted per read, and at the read's end (both mates of a pair) each rule (c, T, m) is applied and tallied, for all
// reads and for those protal writes records of (by the taxon of the best record). Nothing the run writes changes.
//
// The same rules are also tallied for a stride pick of 1/c of a read's cores instead of a hash sample: every c-th core in
// read order (StridePicked), each mate (each chunk of a long read) on its own. A gate looking its sample up in the full
// index (--sparse_first) can pick any cores of a read; a stride spaces them evenly over it, where a hash sample leaves
// gaps by chance (docs/claude/2026-10-13-stride-pick).
//
// Also counted, per mate of a short read: the best score among its first k lookups (in read order, of the cores in the
// index), by whether the mate aligns, its fragment does without it, or neither: whether a mate whose first lookups match
// only by chance could be given up (docs/claude/2026-10-12-seeding-profile, "Giving up on chance matches").
#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include "robin_map.h"

namespace protal::stage1_probe {
    constexpr size_t kLevels = 6;     // c = 1, 2, 4, 8, 16, 32
    constexpr uint32_t kMaxScore = 16;  // a flex cell's 16 flank bases
    constexpr size_t kScores = kMaxScore + 1;
    constexpr std::array<uint8_t, 7> kMinScores{ 10, 11, 12, 13, 14, 15, 16 };  // T
    constexpr std::array<uint8_t, 4> kMinHits{ 1, 2, 3, 5 };                    // m
    constexpr std::array<size_t, 5> kFirst{ 1, 2, 4, 8, SIZE_MAX };             // k (SIZE_MAX: all of the mate's)
    constexpr size_t kTaxonMinHits = 2;  // the taxa table's rules: m = 1 and 2
    constexpr size_t kRules = kLevels * kMinScores.size() * kMinHits.size();
    constexpr size_t kTaxonRules = kLevels * kMinScores.size() * kTaxonMinHits;

    // A mate's state for the first-lookups table: its fragment writes no record; the fragment does, but not with this
    // mate; the mate aligns. A single-end read is the first or the last.
    enum MateState : size_t { kFragmentUnaligned = 0, kOtherMateAligned = 1, kMateAligned = 2, kStates = 3 };

    // The level of a core: the trailing zero bits of its hash (splitmix64's finaliser), at most kLevels - 1. A core is in
    // the sparse index of 1/c, c = 2^j, if its level is j or more, so the sparser indexes are subsets of the denser.
    inline uint8_t Level(uint64_t core) {
        uint64_t z = core + 0x9e3779b97f4a7c15ull;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebull;
        z ^= z >> 31;
        int const zeros = z == 0 ? 64 : __builtin_ctzll(z);
        return static_cast<uint8_t>(std::min<int>(zeros, static_cast<int>(kLevels) - 1));
    }

    // How a rule samples a read's cores: by their hash (Level), or every c-th in read order (StridePicked).
    enum class Pick { hash, stride };

    // The stride pick of 1/c of a read's `cores` cores (c a power of two): every c-th, from the first at StrideOffset, the
    // cores left over shared between the read's two ends (one more at the end where they are odd).
    inline size_t StrideOffset(size_t cores, size_t c) {
        return cores == 0 ? 0 : ((cores - 1) & (c - 1)) / 2;
    }
    inline bool StridePicked(size_t index, size_t cores, size_t c) {
        size_t const offset = StrideOffset(cores, c);
        return index < cores && index >= offset && ((index - offset) & (c - 1)) == 0;
    }

    inline size_t RuleIndex(size_t level, size_t min_score, size_t min_hits) {
        return (level * kMinScores.size() + min_score) * kMinHits.size() + min_hits;
    }

    inline size_t TaxonRuleIndex(size_t level, size_t min_score, size_t min_hits) {
        return (level * kMinScores.size() + min_score) * kTaxonMinHits + min_hits;
    }

    // One read's (a mate's) lookups of cores in the index: counted by level and best score, by stride pick and best
    // score, and in read order.
    struct Mate {
        static constexpr uint32_t kNoIndex = UINT32_MAX;

        std::array<std::array<uint32_t, kScores>, kLevels> counts{};
        // stride_counts[j][s]: lookups of the cores the stride of c = 2^j picks (StridePicked) scoring s.
        std::array<std::array<uint32_t, kScores>, kLevels> stride_counts{};
        std::vector<std::pair<uint32_t, uint8_t>> lookups;  // (read position, best score)
        uint32_t cores = 0;  // the cores of the read (chunk) whose lookups Add counts next (Begin)

        // The read's (a long read's chunk's) cores, those looked up or not, before its lookups are added: the stride's
        // indexes are among these.
        void Begin(uint32_t read_cores) { cores = read_cores; }

        // A lookup: its core's level, its best score, its read position, and its core's index among the read's (chunk's)
        // cores (Begin; kNoIndex: in no stride pick).
        void Add(uint8_t level, uint32_t score, uint32_t read_pos, uint32_t index = kNoIndex) {
            uint8_t const s = static_cast<uint8_t>(std::min(score, kMaxScore));
            counts[std::min<size_t>(level, kLevels - 1)][s]++;
            lookups.emplace_back(read_pos, s);
            if (index == kNoIndex) return;
            for (size_t j = 0; j < kLevels; j++) stride_counts[j][s] += StridePicked(index, cores, size_t{1} << j);
        }

        void Clear() {
            for (auto& row : counts) row.fill(0);
            for (auto& row : stride_counts) row.fill(0);
            lookups.clear();
            cores = 0;
        }
    };

    // The tallies of a thread, joined into Global() at the end of a sample's alignment.
    class Tally {
    public:
        Mate& GetMate(size_t mate) { return m_mates[mate]; }

        // The end of a read: both mates of a pair (`mates` 2), else the first. aligned: protal writes records of it;
        // taxa: the taxa of its best records (one for a short read, one per gene of a long read); mate_aligned: which
        // mates the best record has. short_read: its mates go into the first-lookups table. The mates are cleared.
        void EndRead(size_t mates, bool aligned, std::vector<uint32_t> const& taxa, std::array<bool, 2> mate_aligned, bool short_read) {
            // hits[j][s]: lookups of level j or more scoring s or more, over the read's mates.
            std::array<std::array<uint32_t, kScores + 1>, kLevels + 1> hits{};
            for (size_t j = kLevels; j-- > 0;) {
                for (size_t s = kScores; s-- > 0;) {
                    uint32_t here = 0;
                    for (size_t m = 0; m < mates; m++) here += m_mates[m].counts[j][s];
                    hits[j][s] = here + hits[j + 1][s] + hits[j][s + 1] - hits[j + 1][s + 1];
                }
            }
            // stride_hits[j][s]: lookups of the cores the stride of c = 2^j picks scoring s or more (the strides are not
            // nested as the levels are).
            std::array<std::array<uint32_t, kScores + 1>, kLevels> stride_hits{};
            for (size_t j = 0; j < kLevels; j++) {
                for (size_t s = kScores; s-- > 0;) {
                    uint32_t here = 0;
                    for (size_t m = 0; m < mates; m++) here += m_mates[m].stride_counts[j][s];
                    stride_hits[j][s] = here + stride_hits[j][s + 1];
                }
            }
            m_reads++;
            m_aligned += aligned;
            // The lookups a gate with a rule makes of the read: its sample's, and all of them if it keeps the read.
            uint32_t const all = hits[0][0];
            m_lookups += all;
            for (size_t j = 0; j < kLevels; j++) {
                for (size_t t = 0; t < kMinScores.size(); t++) {
                    uint32_t const n = hits[j][kMinScores[t]];
                    uint32_t const n_stride = stride_hits[j][kMinScores[t]];
                    for (size_t h = 0; h < kMinHits.size(); h++) {
                        size_t const i = RuleIndex(j, t, h);
                        bool const kept = n >= kMinHits[h];
                        m_kept[i] += kept;
                        m_kept_aligned[i] += kept && aligned;
                        m_left[i] += kept ? all : hits[j][0];
                        bool const kept_stride = n_stride >= kMinHits[h];
                        m_kept_stride[i] += kept_stride;
                        m_kept_aligned_stride[i] += kept_stride && aligned;
                        m_left_stride[i] += kept_stride ? all : stride_hits[j][0];
                    }
                }
            }
            if (aligned) {
                for (uint32_t taxid : taxa) {
                    auto& row = m_taxa[taxid];
                    row[0]++;
                    for (size_t j = 0; j < kLevels; j++) {
                        for (size_t t = 0; t < kMinScores.size(); t++) {
                            uint32_t const n = hits[j][kMinScores[t]];
                            uint32_t const n_stride = stride_hits[j][kMinScores[t]];
                            for (size_t h = 0; h < kTaxonMinHits; h++) {
                                row[1 + TaxonRuleIndex(j, t, h)] += n >= kMinHits[h];
                                row[1 + kTaxonRules + TaxonRuleIndex(j, t, h)] += n_stride >= kMinHits[h];
                            }
                        }
                    }
                }
            }
            if (short_read) {
                for (size_t m = 0; m < mates; m++) {
                    auto& lookups = m_mates[m].lookups;
                    std::sort(lookups.begin(), lookups.end());
                    size_t const state = !aligned ? kFragmentUnaligned : mate_aligned[m] ? kMateAligned : kOtherMateAligned;
                    for (size_t f = 0; f < kFirst.size(); f++) {
                        uint8_t best = 0;
                        for (size_t i = 0; i < lookups.size() && i < kFirst[f]; i++) best = std::max(best, lookups[i].second);
                        m_first[f][best][state]++;
                    }
                }
                m_first_mates += mates;
            }
            for (size_t m = 0; m < mates; m++) m_mates[m].Clear();
        }

        void Join(Tally const& other) {
            m_reads += other.m_reads;
            m_aligned += other.m_aligned;
            m_first_mates += other.m_first_mates;
            m_lookups += other.m_lookups;
            for (size_t i = 0; i < kRules; i++) {
                m_kept[i] += other.m_kept[i];
                m_kept_aligned[i] += other.m_kept_aligned[i];
                m_left[i] += other.m_left[i];
                m_kept_stride[i] += other.m_kept_stride[i];
                m_kept_aligned_stride[i] += other.m_kept_aligned_stride[i];
                m_left_stride[i] += other.m_left_stride[i];
            }
            for (auto const& [taxid, row] : other.m_taxa) {
                auto& mine = m_taxa[taxid];
                for (size_t i = 0; i < row.size(); i++) mine[i] += row[i];
            }
            for (size_t f = 0; f < kFirst.size(); f++) {
                for (size_t s = 0; s < kScores; s++) {
                    for (size_t st = 0; st < kStates; st++) m_first[f][s][st] += other.m_first[f][s][st];
                }
            }
        }

        void Clear() { *this = Tally(); }

        uint64_t Reads() const { return m_reads; }
        uint64_t Aligned() const { return m_aligned; }

        // The share of the reads a rule keeps, and of those with records (recall); `pick`: its cores by hash or stride.
        std::pair<double, double> Rule(size_t level, size_t min_score, size_t min_hits, Pick pick = Pick::hash) const {
            size_t const i = RuleIndex(level, min_score, min_hits);
            uint64_t const kept = pick == Pick::hash ? m_kept[i] : m_kept_stride[i];
            uint64_t const kept_aligned = pick == Pick::hash ? m_kept_aligned[i] : m_kept_aligned_stride[i];
            return { m_reads ? static_cast<double>(kept) / static_cast<double>(m_reads) : 0.0,
                     m_aligned ? static_cast<double>(kept_aligned) / static_cast<double>(m_aligned) : 0.0 };
        }

        // The share of the run's lookups (of cores in the index) a gate with a rule makes: every read's sample's, and the
        // rest of the reads it keeps.
        double LookupsLeft(size_t level, size_t min_score, size_t min_hits, Pick pick = Pick::hash) const {
            size_t const i = RuleIndex(level, min_score, min_hits);
            uint64_t const left = pick == Pick::hash ? m_left[i] : m_left_stride[i];
            return m_lookups ? static_cast<double>(left) / static_cast<double>(m_lookups) : 0.0;
        }

        // Every rule of a pick: c, T, m, the reads, those kept, their share, the reads with records, those kept, the
        // recall, and the share of the lookups left (LookupsLeft; since 2026-10-13).
        void WriteRules(std::ostream& os, Pick pick = Pick::hash) const {
            auto const& kept = pick == Pick::hash ? m_kept : m_kept_stride;
            auto const& kept_aligned = pick == Pick::hash ? m_kept_aligned : m_kept_aligned_stride;
            os << "c\tmin_score\tmin_hits\treads\tkept\tshare_kept\taligned\taligned_kept\trecall\tlookups_left\n";
            for (size_t j = 0; j < kLevels; j++) {
                for (size_t t = 0; t < kMinScores.size(); t++) {
                    for (size_t h = 0; h < kMinHits.size(); h++) {
                        size_t const i = RuleIndex(j, t, h);
                        auto const [share, recall] = Rule(j, t, h, pick);
                        os << (1u << j) << '\t' << int(kMinScores[t]) << '\t' << int(kMinHits[h]) << '\t' << m_reads << '\t' << kept[i]
                           << '\t' << Fraction(share) << '\t' << m_aligned << '\t' << kept_aligned[i] << '\t' << Fraction(recall)
                           << '\t' << Fraction(LookupsLeft(j, t, h, pick)) << '\n';
                    }
                }
            }
        }

        // Per taxon of the best records: its reads with records, and those each rule (m = 1 and 2) keeps, the hash
        // sample's (columns c<c>_t<T>_m<m>) then the stride's (s<c>_t<T>_m<m>); `name`, if given, adds the taxon's name.
        // Taxa in taxid order.
        template<typename Name>
        void WriteTaxa(std::ostream& os, Name const& name) const {
            os << "taxid\tname\taligned";
            for (char const pick : { 'c', 's' }) {
                for (size_t j = 0; j < kLevels; j++) {
                    for (size_t t = 0; t < kMinScores.size(); t++) {
                        for (size_t h = 0; h < kTaxonMinHits; h++) os << '\t' << pick << (1u << j) << "_t" << int(kMinScores[t]) << "_m" << int(kMinHits[h]);
                    }
                }
            }
            os << '\n';
            std::map<uint32_t, TaxonRow const*> sorted;
            for (auto const& [taxid, row] : m_taxa) sorted.emplace(taxid, &row);
            for (auto const& [taxid, row] : sorted) {
                os << taxid << '\t' << name(taxid);
                for (uint32_t n : *row) os << '\t' << n;
                os << '\n';
            }
        }

        // Per k and best score of a mate's first k lookups: its mates, by state (MateState).
        void WriteFirstLookups(std::ostream& os) const {
            os << "first_lookups\tbest_score\tmates_fragment_unaligned\tmates_other_mate_aligned\tmates_aligned\n";
            for (size_t f = 0; f < kFirst.size(); f++) {
                for (size_t s = 0; s < kScores; s++) {
                    auto const& row = m_first[f][s];
                    if (row[0] + row[1] + row[2] == 0) continue;
                    os << (kFirst[f] == SIZE_MAX ? std::string("all") : std::to_string(kFirst[f])) << '\t' << s << '\t' << row[0]
                       << '\t' << row[1] << '\t' << row[2] << '\n';
                }
            }
        }

        // Of the mates (short reads), the share whose best score among their first k lookups is at most `score`, of all
        // mates and of those that align (and of the mates whose fragment aligns).
        std::array<double, 3> GiveUp(size_t first, uint8_t score) const {
            std::array<uint64_t, kStates> low{}, all{};
            for (size_t s = 0; s < kScores; s++) {
                for (size_t st = 0; st < kStates; st++) {
                    all[st] += m_first[first][s][st];
                    if (s <= score) low[st] += m_first[first][s][st];
                }
            }
            auto share = [](uint64_t a, uint64_t b) { return b ? static_cast<double>(a) / static_cast<double>(b) : 0.0; };
            return { share(low[0] + low[1] + low[2], all[0] + all[1] + all[2]), share(low[2], all[2]),
                     share(low[1] + low[2], all[1] + all[2]) };
        }

        uint64_t FirstLookupMates() const { return m_first_mates; }

    private:
        static std::string Fraction(double value) {
            std::ostringstream os;
            os << std::fixed << std::setprecision(6) << value;
            return os.str();
        }

        // A taxon's reads with records, then those the hash sample's rules keep, then the stride's.
        using TaxonRow = std::array<uint32_t, 1 + 2 * kTaxonRules>;

        std::array<Mate, 2> m_mates;
        uint64_t m_reads = 0, m_aligned = 0, m_first_mates = 0, m_lookups = 0;
        std::array<uint64_t, kRules> m_kept{}, m_kept_aligned{}, m_left{}, m_kept_stride{}, m_kept_aligned_stride{}, m_left_stride{};
        tsl::robin_map<uint32_t, TaxonRow> m_taxa;
        std::array<std::array<std::array<uint64_t, kStates>, kScores>, kFirst.size()> m_first{};
    };

    // The sample's tallies, joined from the threads' (under the loops' critical section), written and cleared by
    // RunProtal after each sample.
    inline Tally& Global() {
        static Tally tally;
        return tally;
    }

    // Whether the run probes (--stage1_probe): set by RunProtal before the samples are aligned, read by Classify's loops.
    inline bool& EnabledFlag() {
        static bool enabled = false;
        return enabled;
    }
    inline bool Enabled() { return EnabledFlag(); }

    // "c=8 T=14 m=1: 3.10% kept, recall 99.60%" (a stride pick's "stride c=8 ..."): the log's rules.
    inline std::string RuleText(Tally const& tally, size_t level, size_t min_score, size_t min_hits, Pick pick = Pick::hash) {
        auto const [share, recall] = tally.Rule(level, min_score, min_hits, pick);
        std::ostringstream os;
        os << std::fixed << std::setprecision(2) << (pick == Pick::stride ? "stride c=" : "c=") << (1u << level) << " T=" << int(kMinScores[min_score]) << " m="
           << int(kMinHits[min_hits]) << ": " << 100.0 * share << "% kept, recall " << 100.0 * recall << '%';
        return os.str();
    }
}
