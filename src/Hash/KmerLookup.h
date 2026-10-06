//
// Created by fritsche on 18/08/22.
//

#pragma once

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <algorithm>
#include <tuple>
#include "Seedmap.h"
#include "FlexScan.h"
#include "Constants.h"
#include "TargetClones.h"

namespace protal {
    // KmerPutter is a Metatemplate programming "Interface"  and needs to satisfy the following functions.
    // KmerPutter needs to handle OMP collisions or have the map datastructure deal with it internally
    // -> Put function:  Put(size_t key, size_t taxid, size_t geneid, size_t genepos)
    // -> Boolean that says if it needs pre_processing
    //    template <typename T, typename S=size_t>
    //    concept KmerPutter = requires (T t, S s) {
    //        { t.Put(s, s, s, s) } -> std::same_as<void>;
    //        { t.FirstPut(s, s, s, s) } -> std::same_as<void>;
    //    };

    struct LookupResult {
        uint32_t taxid = UINT32_MAX;
        uint32_t geneid = UINT32_MAX;
        uint32_t genepos = UINT32_MAX;
        uint16_t readpos = UINT16_MAX;
        bool unique = false;
        bool unique_dist_two = false;

        LookupResult() {};

        LookupResult(uint32_t taxid, uint32_t geneid, uint32_t genepos, uint32_t readpos, bool unique=false, bool unique_dist_two=false) :
                taxid(taxid), geneid(geneid), genepos(genepos), readpos(readpos), unique(unique), unique_dist_two(unique_dist_two) {};

        inline bool SameTaxon(LookupResult const& other) const {
            return taxid == other.taxid;
        }

        inline bool SameGene(LookupResult const& other) const {
            return geneid == other.geneid;
        }

        inline bool FromSameSequence(LookupResult const& other) const {
            return geneid == other.geneid && taxid == other.taxid;
        }

        inline bool Equals(LookupResult const& other) const {
            return geneid == other.geneid && taxid == other.taxid && genepos == other.genepos;
        }

        static inline int IntComp(uint32_t const& a, uint32_t const& b) {
            return (a < b) * 1 + (a > b) * -1;
        }

//        inline int operator <=> (LookupResult const& other) const {
//            bool taxid_comp = IntComp(this->taxid, other.taxid);
//            if (taxid_comp != 0) {
//                return taxid_comp;
//            }
//            return IntComp(this->geneid, other.geneid);
//        }

        inline bool operator > (LookupResult const& other) const {
            if (this->taxid != other.taxid) {
                return this->taxid > other.taxid;
            }
            if (this->geneid != other.geneid) {
                return this->geneid > other.geneid;
            }
            return false;
        };

        template<uint64_t taxshift, uint64_t geneshift>
        inline uint64_t ToUINT64_t() const {
            return (static_cast<uint64_t>(taxid) << taxshift) | (static_cast<uint64_t>(geneid) << geneshift) | static_cast<uint64_t>(readpos);
        }

        inline uint64_t ToUINT64_t2() const {
            return (static_cast<uint64_t>(taxid) << 40llu) | (static_cast<uint64_t>(geneid) << 20llu) | static_cast<uint64_t>(readpos);
        }

        inline bool operator < (LookupResult const& other) const {
            if (this->taxid != other.taxid) {
                return this->taxid < other.taxid;
            }
            if (this->geneid != other.geneid) {
                return this->geneid < other.geneid;
            }
            return false;
        };
//        inline bool operator > (LookupResult const& other) const = default;
        inline int operator == (LookupResult const& other) const {
            return this->taxid == other.taxid && this->geneid == other.geneid;
        }

        static bool SortByReadComparator(LookupResult const& a, LookupResult const& b) {
            if (a.taxid != b.taxid) {
                return a.taxid < b.taxid;
            }
            if (a.geneid != b.geneid) {
                return a.geneid < b.geneid;
            }
            if (a.readpos != b.readpos) {
                return a.readpos < b.readpos;
            }
            return false;
        }

        static bool SortByReadComparator2(LookupResult const& a, LookupResult const& b) {
            return a.ToUINT64_t<40,20>() < b.ToUINT64_t<40,20>();
        }

        // A read's seeds are sorted (SharedSeeds::Sort, ChainAnchorFinder.h) by taxon, gene and read position, as SortByReadComparator
        // orders them, and then by gene position and the two flags: every field takes part, so any sort gives the same
        // order. (std::sort with SortByReadComparator left seeds equal in the first three, a k-mer's hits at several
        // places of one gene, in an order of the standard library's making, and the first of a gene's seeds starts its
        // anchor.) All fields packed into one 128-bit key, which also sorts in half the time of the comparator.
        __extension__ using SortKeyType = unsigned __int128;
        SortKeyType SortKey() const {
            return (SortKeyType{ taxid } << 82) | (SortKeyType{ geneid } << 50) | (SortKeyType{ readpos } << 34) |
                   (SortKeyType{ genepos } << 2) | (SortKeyType{ unique } << 1) | SortKeyType{ unique_dist_two };
        }
        static LookupResult FromSortKey(SortKeyType key) {
            LookupResult seed;
            seed.taxid = static_cast<uint32_t>(key >> 82);
            seed.geneid = static_cast<uint32_t>(key >> 50);
            seed.readpos = static_cast<uint16_t>(key >> 34);
            seed.genepos = static_cast<uint32_t>(key >> 2);
            seed.unique = (key >> 1) & 1;
            seed.unique_dist_two = key & 1;
            return seed;
        }

        std::string ToString() const {
            std::string str;
            str += "[LUR: ";
            str += std::to_string(taxid) + ", ";
            str += std::to_string(geneid) + ", RPOS:";
            str += std::to_string(genepos) + ", QPOS:";
            str += std::to_string(readpos) + "]";
            return str;
        }
    };

    // A k-mer's values in the index (Seedmap::PackedBlock: where its entries and flex cells lie and how
    // many), with the k-mer's own flex part and its position in the read.
    struct LookupPointer : Seedmap::PackedBlock {
        uint32_t flex_key = 0;
        uint32_t read_pos = 0;
        bool retrieved = false;

        std::string ToString() {
            std::string str;
            str += "[" + std::to_string(read_pos) + ", " + std::to_string(size) + ", " + KmerUtils::ToString(flex_key, 32) + "]";
            return str;
        }
    };



    class KmerLookupSM {
        std::vector<LookupResult> m_results;
        std::vector<LookupPointer> m_lookups;
//        std::shared_ptr<Seedmap> m_sm;
        Seedmap& m_sm;

        std::vector<uint8_t> flex_vector;  // the last block's scores (flex_scan::BestAvx2, ScoreScalar), grown as needed
        std::vector<uint32_t> m_tie_masks;  // its best cells, 32 to a word (flex_scan::TiesAvx2)

        size_t m_flex_k = 16;
        size_t m_flex_k_half = m_flex_k/2;
        size_t m_flex_k_bits = m_flex_k*2;
        size_t m_max_ubiquity = 16;

        size_t m_taxid{};
        size_t m_geneid{};
        size_t m_genepos{};
        bool m_unique{};
        bool m_unique_dist_two{};

        LookupPointer m_lookup_tmp;
    public:
//        inline void Load(std::istream& ifs) {
////            m_sm = std::make_shared<Seedmap>(Seedmap());
////            m_sm->Load(ifs);
//        }
        using RecoverySet = tsl::robin_set<uint64_t>;

        KmerLookupSM(Seedmap& map, size_t max_ubiquity) :
                m_sm(map), m_flex_k(map.m_flex_k), m_flex_k_bits(map.m_flex_k_bits), m_flex_k_half(map.m_flex_k/2),
                m_max_ubiquity(max_ubiquity) {
        };

        KmerLookupSM(Seedmap& map) :
            m_sm(map), m_flex_k(map.m_flex_k), m_flex_k_bits(map.m_flex_k_bits), m_flex_k_half(map.m_flex_k/2),
            m_max_ubiquity(map.max_key_multiplicity) {
        };

        KmerLookupSM(KmerLookupSM const& other) :
                m_sm(other.m_sm), m_flex_k(other.m_sm.m_flex_k), m_flex_k_bits(other.m_sm.m_flex_k_bits),
                m_max_ubiquity(other.m_max_ubiquity), m_flex_k_half(other.m_sm.m_flex_k/2) {
        };

        void Clear() {
            m_lookups.clear();
        }

        // What Get and GetFromLookup will read, fetched ahead (Seedmap::PrefetchKey): the key map's block of a
        // k-mer, and the first cache line of a lookup's values and of its flex cells.
        inline void PrefetchKey(size_t kmer) {
            m_sm.PrefetchKey(kmer);
        }

        // What GetFromLookup reads first: the flex cells (if the key has them) and the entries.
        static inline void PrefetchValues(LookupPointer const& pointers) {
            if (pointers.flex != nullptr) __builtin_prefetch(pointers.flex);
            __builtin_prefetch(pointers.entries);
        }

        // Entry i of a lookup's (or a core's) values as a ValueEntry (Seedmap::EntryValue).
        inline ValueEntry Entry(Seedmap::PackedBlock const& pointers, uint32_t i) const {
            ValueEntry entry;
            entry.value = m_sm.EntryValue(pointers, i);
            return entry;
        }

        inline void Get(std::vector<LookupPointer>& result, size_t &kmer, uint32_t readpos) {
            if (!m_sm.GetPacked(kmer, m_lookup_tmp)) return;
            m_lookup_tmp.flex_key = static_cast<uint32_t>(m_sm.FlexKey(kmer));
            m_lookup_tmp.read_pos = readpos;
            result.emplace_back(m_lookup_tmp);
        }

        // The seeds of a lookup: the entries whose flex cells share the most bases with the read's (all entries of a key
        // with one value, which has no flex cells), in entry order, flagged unique only where the whole k-mer matches.
        // None if more than m_max_ubiquity cells share the best score; returns false for such a lookup (too ubiquitous).
        PROTAL_CLONE_V3 inline bool GetFromLookup(LookupList& result, LookupPointer& pointers) {
            if (pointers.flex == nullptr) {
                for (uint32_t i = 0; i < pointers.size; i++) {
                    Entry(pointers, i).Get(m_taxid, m_geneid, m_genepos, m_unique, m_unique_dist_two);
                    result.emplace_back( m_taxid, m_geneid, m_genepos, pointers.read_pos + m_flex_k_half, false, false );
                }
                return true;
            }
            // The cells with the best score as bit masks (m_tie_masks, 32 cells to a word), then their entries in order.
            uint32_t const size = pointers.size;
            uint32_t const words = static_cast<uint32_t>(flex_scan::TieWords(size));
            if (m_tie_masks.size() < words) m_tie_masks.resize(words);
            uint32_t max = 0, max_count = 0;
            if (flex_scan::Avx2Enabled().load(std::memory_order_relaxed)) {
                // Every cell's score and the best in one pass, the masks and their count in another (FlexScan.h).
                size_t const bytes = std::max<size_t>(size + flex_scan::kScorePadding, flex_scan::TieScoreBytes(size));
                if (flex_vector.size() < bytes) flex_vector.resize(bytes);
                max = flex_scan::BestAvx2(pointers, pointers.flex_key, flex_vector.data());
                max_count = flex_scan::TiesAvx2(flex_vector.data(), size, max, m_tie_masks.data());
                if (max_count > m_max_ubiquity) return false;
            } else {
                // Every cell's score, the best and how many have it, then the masks one cell at a time.
                if (flex_vector.size() < size) flex_vector.resize(size);
                std::tie(max, max_count) = flex_scan::ScoreScalar(pointers, pointers.flex_key, flex_vector.data());
                if (max_count > m_max_ubiquity) return false;
                std::fill_n(m_tie_masks.begin(), words, 0u);
                for (uint32_t i = 0; i < size; i++) m_tie_masks[i / 32] |= static_cast<uint32_t>(flex_vector[i] == max) << (i % 32);
            }
            bool const exact = max == m_sm.m_flex_k;
            for (uint32_t w = 0; w < words; w++) {
                for (uint32_t m = m_tie_masks[w]; m != 0; m &= m - 1) {
                    ValueEntry const entry = Entry(pointers, 32 * w + static_cast<uint32_t>(__builtin_ctz(m)));
                    entry.Get(m_taxid, m_geneid, m_genepos, m_unique, m_unique_dist_two);
                    m_unique &= exact;
                    m_unique_dist_two &= exact;

                    if (m_taxid == 0) {
                        // Debug
                        std::cout << ValueEntry(entry).ToString() << std::endl;
                        for (auto& e : result) {
                            std::cout << e.ToString() << std::endl;
                        }
                        continue;
                    }

                    result.emplace_back( m_taxid, m_geneid, m_genepos, pointers.read_pos + m_flex_k_half, m_unique, m_unique_dist_two );
                }
            }
            return true;
        }

        inline void RecoverFromLookup(LookupList& result, LookupPointer& pointers, RecoverySet& recovery) {
            for (uint32_t i = 0; i < pointers.size; i++) {
                ValueEntry entry = Entry(pointers, i);
                if (recovery.contains(entry.MaskPosition())) {
                    auto [taxid, geneid, genepos] = entry.Get();
                    result.emplace_back(LookupResult(taxid, geneid, genepos, pointers.read_pos + m_flex_k_half));
                    recovery.erase(recovery.find(entry.MaskPosition()));
                    return;
                }
            }
        }

        // Whether a k-mer core occurs only once in the index; its values in `block` (one entry, no flex
        // cells: see m_flex_threshold, so GetExact returns nothing for it).
        inline bool GetSingleEntry(size_t &kmer, Seedmap::PackedBlock& block) const {
            return m_sm.GetPacked(kmer, block) && block.flex == nullptr && block.size == 1;
        }

        // The entries (their indices in `block`, the core's values) of a core with several values whose
        // whole k-mer is `kmer` (equal flex parts), none if more than m_max_ubiquity are; appended to
        // `entries`, which must be empty. For --build's uniqueness check, which only acts on whole k-mers.
        // Returns the flex parts compared.
        inline size_t GetExact(size_t &kmer, Seedmap::PackedBlock& block, std::vector<uint32_t>& entries) const {
            if (!m_sm.GetPacked(kmer, block) || block.flex == nullptr) return 0;
            uint32_t const flex_key = static_cast<uint32_t>(m_sm.FlexKey(kmer));
            for (uint32_t i = 0; i < block.size; i++) {
                if (Seedmap::FlexCell(block, i) == flex_key) entries.emplace_back(i);
            }
            if (entries.size() > m_max_ubiquity) entries.clear();
            return block.size;
        }

        // Get and GetFromLookup in one.
        inline void Get(LookupList& result, size_t &kmer, uint32_t readpos, RecoverySet* choose=nullptr, LookupList* recovered_results=nullptr) {
            (void)choose;
            (void)recovered_results;
            if (!m_sm.GetPacked(kmer, m_lookup_tmp)) return;
            m_lookup_tmp.flex_key = static_cast<uint32_t>(m_sm.FlexKey(kmer));
            m_lookup_tmp.read_pos = readpos;
            GetFromLookup(result, m_lookup_tmp);
        }
    };
}
