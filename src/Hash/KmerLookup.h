//
// Created by fritsche on 18/08/22.
//

#pragma once

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include "Seedmap.h"
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

        std::vector<uint16_t> flex_vector;

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

        PROTAL_CLONE_V3 inline void GetFromLookup(LookupList& result, LookupPointer& pointers) {
            if (pointers.flex != nullptr) {
                flex_vector.clear();
                auto max = 0;
                auto max_count = 0;
                for (uint32_t i = 0; i < pointers.size; i++) {
                    auto sim = Seedmap::Similarity(Seedmap::FlexCell(pointers, i), pointers.flex_key);
                    flex_vector.emplace_back(sim);
                    if (sim > max)  {
                        max = sim;
                        max_count = 0;
                    }
                    max_count += (sim == max);
                }

                if (max_count > m_max_ubiquity) {
                    return;
                }

                for (uint32_t i = 0; i < flex_vector.size(); i++) {
                    if (flex_vector[i] == max) {
                        ValueEntry const entry = Entry(pointers, i);
                        entry.Get(m_taxid, m_geneid, m_genepos, m_unique, m_unique_dist_two);
                        m_unique &= max == m_sm.m_flex_k;
                        m_unique_dist_two &= max == m_sm.m_flex_k;

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
            } else {
                for (uint32_t i = 0; i < pointers.size; i++) {
                    Entry(pointers, i).Get(m_taxid, m_geneid, m_genepos, m_unique, m_unique_dist_two);
                    result.emplace_back( m_taxid, m_geneid, m_genepos, pointers.read_pos + m_flex_k_half, false, false );
                }
            }
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
