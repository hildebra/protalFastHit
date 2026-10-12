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

    inline size_t RuleIndex(size_t level, size_t min_score, size_t min_hits) {
        return (level * kMinScores.size() + min_score) * kMinHits.size() + min_hits;
    }

    inline size_t TaxonRuleIndex(size_t level, size_t min_score, size_t min_hits) {
        return (level * kMinScores.size() + min_score) * kTaxonMinHits + min_hits;
    }

    // One read's (a mate's) lookups of cores in the index: counted by level and best score, and in read order.
    struct Mate {
        std::array<std::array<uint32_t, kScores>, kLevels> counts{};
        std::vector<std::pair<uint32_t, uint8_t>> lookups;  // (read position, best score)

        void Add(uint8_t level, uint32_t score, uint32_t read_pos) {
            uint8_t const s = static_cast<uint8_t>(std::min(score, kMaxScore));
            counts[std::min<size_t>(level, kLevels - 1)][s]++;
            lookups.emplace_back(read_pos, s);
        }

        void Clear() {
            for (auto& row : counts) row.fill(0);
            lookups.clear();
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
            m_reads++;
            m_aligned += aligned;
            for (size_t j = 0; j < kLevels; j++) {
                for (size_t t = 0; t < kMinScores.size(); t++) {
                    uint32_t const n = hits[j][kMinScores[t]];
                    for (size_t h = 0; h < kMinHits.size(); h++) {
                        bool const kept = n >= kMinHits[h];
                        m_kept[RuleIndex(j, t, h)] += kept;
                        m_kept_aligned[RuleIndex(j, t, h)] += kept && aligned;
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
                            for (size_t h = 0; h < kTaxonMinHits; h++) row[1 + TaxonRuleIndex(j, t, h)] += n >= kMinHits[h];
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
            for (size_t i = 0; i < kRules; i++) {
                m_kept[i] += other.m_kept[i];
                m_kept_aligned[i] += other.m_kept_aligned[i];
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

        // The share of the reads a rule keeps, and of those with records (recall).
        std::pair<double, double> Rule(size_t level, size_t min_score, size_t min_hits) const {
            size_t const i = RuleIndex(level, min_score, min_hits);
            return { m_reads ? static_cast<double>(m_kept[i]) / static_cast<double>(m_reads) : 0.0,
                     m_aligned ? static_cast<double>(m_kept_aligned[i]) / static_cast<double>(m_aligned) : 0.0 };
        }

        // Every rule: c, T, m, the reads, those kept, their share, the reads with records, those kept, the recall.
        void WriteRules(std::ostream& os) const {
            os << "c\tmin_score\tmin_hits\treads\tkept\tshare_kept\taligned\taligned_kept\trecall\n";
            for (size_t j = 0; j < kLevels; j++) {
                for (size_t t = 0; t < kMinScores.size(); t++) {
                    for (size_t h = 0; h < kMinHits.size(); h++) {
                        size_t const i = RuleIndex(j, t, h);
                        auto const [share, recall] = Rule(j, t, h);
                        os << (1u << j) << '\t' << int(kMinScores[t]) << '\t' << int(kMinHits[h]) << '\t' << m_reads << '\t' << m_kept[i]
                           << '\t' << Fraction(share) << '\t' << m_aligned << '\t' << m_kept_aligned[i] << '\t' << Fraction(recall) << '\n';
                    }
                }
            }
        }

        // Per taxon of the best records: its reads with records, and those each rule (m = 1 and 2) keeps; `name`, if
        // given, adds the taxon's name. Taxa in taxid order.
        template<typename Name>
        void WriteTaxa(std::ostream& os, Name const& name) const {
            os << "taxid\tname\taligned";
            for (size_t j = 0; j < kLevels; j++) {
                for (size_t t = 0; t < kMinScores.size(); t++) {
                    for (size_t h = 0; h < kTaxonMinHits; h++) os << "\tc" << (1u << j) << "_t" << int(kMinScores[t]) << "_m" << int(kMinHits[h]);
                }
            }
            os << '\n';
            std::map<uint32_t, std::array<uint32_t, 1 + kTaxonRules> const*> sorted;
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

        std::array<Mate, 2> m_mates;
        uint64_t m_reads = 0, m_aligned = 0, m_first_mates = 0;
        std::array<uint64_t, kRules> m_kept{}, m_kept_aligned{};
        tsl::robin_map<uint32_t, std::array<uint32_t, 1 + kTaxonRules>> m_taxa;  // [0]: reads with records, then the rules
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

    // "c=8 T=14 m=1: kept 3.1% of the reads, 99.6% of those with records": the log's rules.
    inline std::string RuleText(Tally const& tally, size_t level, size_t min_score, size_t min_hits) {
        auto const [share, recall] = tally.Rule(level, min_score, min_hits);
        std::ostringstream os;
        os << std::fixed << std::setprecision(2) << "c=" << (1u << level) << " T=" << int(kMinScores[min_score]) << " m="
           << int(kMinHits[min_hits]) << ": " << 100.0 * share << "% kept, recall " << 100.0 * recall << '%';
        return os.str();
    }
}
