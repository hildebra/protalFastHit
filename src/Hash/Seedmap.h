//
// Created by fritsche on 17/02/2022.
//

#pragma once

#include <Constants.h>
#include <assert.h>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <sys/mman.h>
#include <unistd.h>
#include <bitset>
#include <iostream>
#include <tuple>
#include <fstream>
#include "Utilities.h"
#include "Zstd.h"
#include "IndexCodec.h"
#include "KmerUtils.h"
#include "ReferenceFingerprint.h"
#include "sparse_map.h"
#include <atomic>
#include <bit>
#include <bits/stdc++.h>
#include "protal_config.h"

namespace protal {
    template<uint64_t taxid_bits, uint64_t geneid_bits, uint64_t genepos_bits>
    struct Entry {
        uint64_t value = 0;

        inline void Put(uint64_t taxid, uint64_t geneid, uint64_t genepos) {
            value = 0;
            value |= taxid << (geneid_bits + genepos_bits);
            value |= geneid << genepos_bits;
            value |= genepos;

            // New Flag defaults to one
            value |= 1llu << (taxid_bits + geneid_bits + genepos_bits);
        }

        inline uint64_t MaskPosition() {
            return value & (((1llu << (taxid_bits + geneid_bits)) - 1) << genepos_bits);
        }

        void Get(uint64_t &taxid, uint64_t &geneid, uint64_t &genepos) const {
            taxid = (value >> (geneid_bits + genepos_bits)) & ((1u << taxid_bits) - 1);
            geneid = (value >> (genepos_bits)) & ((1u << geneid_bits) - 1);
            genepos = value & ((1u << genepos_bits) - 1);

            // New flag is one by default
            bool flag1 = (value >> (geneid_bits + genepos_bits + taxid_bits)) & 1;
        }

        void Get(uint64_t &taxid, uint64_t &geneid, uint64_t &genepos, bool& unique, bool& unique_min_distance_two) const {
            taxid = (value >> (geneid_bits + genepos_bits)) & ((1u << taxid_bits) - 1);
            geneid = (value >> (genepos_bits)) & ((1u << geneid_bits) - 1);
            genepos = value & ((1u << genepos_bits) - 1);

            // New flag is one by default
            unique = IsFlagUnique();
            unique_min_distance_two = IsFlagUniqueDistanceMinTwo();
        }

        // Atomic: --build's uniqueness check clears flags from all its threads (with a lock around
        // each, 8 threads were slower than 4). Clearing a bit commutes, so the result does not depend on
        // the order.
        void SetFlagNonUnique() {
            std::atomic_ref<uint64_t>(value).fetch_and(~(1llu << (taxid_bits + geneid_bits + genepos_bits)),
                                                       std::memory_order_relaxed);
        }

        void SetFlagUniqueDistanceMinTwo() {
            value |= (1llu << (taxid_bits + geneid_bits + genepos_bits + 1));
        }

        size_t IsFlagUnique() const {
            return (value >> (geneid_bits + genepos_bits + taxid_bits)) & 1;
        }

        size_t IsFlagUniqueDistanceMinTwo() const {
            return (value >> (geneid_bits + genepos_bits + taxid_bits + 1)) & 1;
        }

        inline bool operator > (Entry const& other) const {
            auto [this_taxid, this_geneid, this_genepos] = Get();
            auto [other_taxid, other_geneid, other_genepos] = other.Get();
            if (this_taxid != other_taxid) {
                return this_taxid > other_taxid;
            }
            if (this_geneid != other_geneid) {
                return this_geneid > other_geneid;
            }
            return false;
        }

        inline bool operator < (Entry const& other) const {
            auto [this_taxid, this_geneid, this_genepos] = Get();
            auto [other_taxid, other_geneid, other_genepos] = other.Get();
            if (this_taxid != other_taxid) {
                return this_taxid < other_taxid;
            }
            if (this_geneid != other_geneid) {
                return this_geneid < other_geneid;
            }
            return false;
        };

//        inline bool operator > (LookupResult const& other) const = default;
        inline int operator == (Entry const& other) const {
            auto [this_taxid, this_geneid, this_genepos] = Get();
            auto [other_taxid, other_geneid, other_genepos] = other.Get();
            return this_taxid == other_taxid && this_geneid == other_geneid;
        }

        std::tuple<uint64_t, uint64_t, uint64_t> Get() const {
            auto taxid = (value >> (geneid_bits + genepos_bits)) & ((1u << taxid_bits) - 1);
            auto geneid = (value >> genepos_bits) & ((1u << geneid_bits) - 1);
            auto genepos = value & ((1u << genepos_bits) - 1);
            return { taxid, geneid, genepos };
        }

        std::string ToString() {
            std::string result;
            auto [taxid, geneid, genepos] = Get();
            result.reserve(60);
            result += "taxid:\t" + std::to_string(taxid) + "\t\t";
            result += "geneid:\t" + std::to_string(geneid) + "\t\t";
            result += "genepos:\t" + std::to_string(genepos) + "\t\t";
            result += "unique: " + std::to_string(IsFlagUnique());
            return result;
        }

        bool Empty() {
            return value == 0;
        }
    };

    using ValueEntry = Entry<SEEDMAP_TAXID_BITS, SEEDMAP_GENEID_BITS, SEEDMAP_GENE_POS_BITS>;
    static_assert(std::is_trivially_copyable_v<ValueEntry>, "index values are read as raw bytes");

    // Scatters an index file's bytes (from zstd::ParallelRead) into the key map and the values;
    // the header before them has already been read. Bytes past the values are dropped (the
    // caller compares the size read with the size expected).
    class IndexSink : public zstd::Sink {
    public:
        IndexSink(uint64_t header_bytes, char* keymap, uint64_t keymap_bytes, char* values, uint64_t values_bytes)
                : m_keymap_begin(header_bytes), m_values_begin(header_bytes + keymap_bytes),
                  m_end(header_bytes + keymap_bytes + values_bytes), m_keymap(keymap), m_values(values) {}

        char* Direct(uint64_t offset, size_t size) override {
            if (offset >= m_keymap_begin && offset + size <= m_values_begin) return m_keymap + (offset - m_keymap_begin);
            if (offset >= m_values_begin && offset + size <= m_end) return m_values + (offset - m_values_begin);
            return nullptr;
        }

        void Copy(uint64_t offset, char const* data, size_t size) override {
            CopyPart(offset, data, size, m_keymap_begin, m_values_begin, m_keymap);
            CopyPart(offset, data, size, m_values_begin, m_end, m_values);
        }

    private:
        static void CopyPart(uint64_t offset, char const* data, size_t size, uint64_t begin, uint64_t end, char* dst) {
            uint64_t const from = std::max(offset, begin), to = std::min(offset + size, end);
            if (from < to) std::memcpy(dst + (from - begin), data + (from - offset), static_cast<size_t>(to - from));
        }

        uint64_t m_keymap_begin, m_values_begin, m_end;
        char* m_keymap;
        char* m_values;
    };

    class SeedmapUtils {
    public:
        template<size_t bits>
        static std::string BitString(uint64_t key) {
            std::bitset<bits> x(key);
            return x.to_string();
        }

        static inline void ParseHeader(std::string &header, size_t &internal_taxid, size_t &marker_gene_id) {
            auto tokens = Utils::split(header.substr(1, std::string::npos), "_");
            assert(tokens.size() > 0);
            internal_taxid = stoull(tokens[0]);
            marker_gene_id = (tokens.size() > 1) * stoull(tokens[1]);
        }
    };

//    │  values table                 │
//    │                               │
//    │                               │
//    ├───────────────┬───────────────┤ ─┐
//    │ 1  TCACACGTC  │ 2 GTCACACGC   │  │
//    ├───────────────┼───────────────┤  │   Here are the flexi-k values stored
//    │ 3  ATGCATGCT  │ 4 CGACTCGGC   │  │   Use a loop and similarity check
//    ├───────────────┼───────────────┤  │   Maybe sth with popcount to find
//    │ 5  ...        │ 6  ...        │  │   best matching k-mer in the value-section
//    ├───────────────┼───────────────┤  │
//    │ 7  ...        │ 8  ...        │  │
// ┌─ ├───────────────┴───────────────┤ ─┘
// │  │ 1 ( taxid, geneid, genepos )  │
// │  ├───────────────────────────────┤
// │  │ 2                             │
// │  ├───────────────────────────────┤
// │  │               .               │
// │  │               .               │
// │  │               .               │
// │  ├───────────────────────────────┤
// │  │ 7                             │
// │  ├───────────────────────────────┤
// │  │ 8                             │
// └─ ├───────────────────────────────┤
//    │                               │


    class Seedmap {
    public:
        // Format 1 ("PRXSEEDM"): magic and version. Format 2 ("PRXSEED2") adds a feature word after
        // the version; protal builds before format 2 reject it as an unknown header.
        static constexpr uint64_t kFileMagic = 0x505258534545444dllu;
        static constexpr uint64_t kFileMagicV2 = 0x5052585345454432llu;
        // Features of an index. A format 1 index has none of them.
        static constexpr uint32_t kFeatureFullSyncmerMask = 1u << 0;       // syncmers compare whole s-mers
        static constexpr uint32_t kFeatureCorrectUniqueTwoFlags = 1u << 1; // "unique at distance >= 2" per entry
        static constexpr uint32_t kFeatureReferenceFingerprint = 1u << 2;  // header carries a ReferenceFingerprint
        static constexpr uint32_t kFeatureCheckedSingleEntries = 1u << 3;  // single-entry k-mers checked for uniqueness
        static constexpr uint32_t kKnownFeatures = kFeatureFullSyncmerMask | kFeatureCorrectUniqueTwoFlags |
                                                   kFeatureReferenceFingerprint | kFeatureCheckedSingleEntries;
        static constexpr uint32_t kIndexVersionMajor = protal_VERSION_MAJOR;
        static constexpr uint32_t kIndexVersionMinor = protal_VERSION_MINOR;
        static constexpr uint32_t kIndexVersionPatch = protal_VERSION_PATCH;
        static constexpr uint32_t kOldestCompatibleIndexVersionMajor = 0;
        static constexpr uint32_t kOldestCompatibleIndexVersionMinor = 3;
        static constexpr uint32_t kOldestCompatibleIndexVersionPatch = 0;

        // If a 15-mer has more than 8 locations, use flexi-k approach
        size_t m_exact_k = 15;
        size_t m_main_bits = m_exact_k * 2;
        size_t m_flex_threshold = 2;
        size_t m_flex_k = 16;
        size_t m_flex_k_half = m_flex_k/2;
        size_t m_flex_k_bits = m_flex_k*2;
        size_t m_main_key_mask = ((1llu << m_main_bits) - 1) << m_flex_k; // First mask, then shift
        size_t m_flex_key_mask_left = ((1llu << m_flex_k) - 1) << (m_flex_k + m_main_bits);
        size_t m_flex_key_mask_right = (1llu << m_flex_k) - 1;

        size_t max_key_multiplicity = 2048; //  This is for key, arbitrary value

    private:

//        using KeyMap_t = uint8_t;
        using KeyMap_t = uint16_t; // CHANGE THIS IN LATEST VERSION
        KeyMap_t* m_keymap = nullptr;
        uint64_t m_keymask = 0b0000000000000000000000000000000000111111111111111111111111111111;
        size_t keymap_size = 1u << m_main_bits;
        size_t keymap_max = m_keymask;

        // This variable defines how many keys are managed by one control block
        // Must be power of two
//        size_t m_keys_per_ctrl_block = 16; // Version 1
        size_t m_keys_per_ctrl_block = 8; // Version 2
//        size_t m_keys_per_ctrl_block = 32; // Flex 16-bit Cells

        // bitshift to find which control block a key is in.
        size_t ctrl_block_frequency_bitshift = log2(m_keys_per_ctrl_block);

        uint64_t ctrl_block_key_mask = (1 << ctrl_block_frequency_bitshift) - 1;

        // (1 << 8) = 256. This is how many fields an 8-bit number can index.
        // Spread the indexing power equally across the number of keys per block
        size_t max_key_ubiquity = (1 << (sizeof(KeyMap_t)*8)) / m_keys_per_ctrl_block;

//        size_t ctrl_block_frequency_bitshift = 4;

        // Bytespace a controlblock takes up
        size_t ctrl_block_byte_size = 8; // Byte
        size_t ctrl_block_cell_size = ctrl_block_byte_size/sizeof(KeyMap_t); // Byte
        size_t ctrl_block_byte_size_shift = log2(ctrl_block_cell_size);
//        size_t ctrl_block_byte_size_shift = 2;

        size_t keymap_size_total = keymap_size + (((keymap_size / m_keys_per_ctrl_block) + 1) * ctrl_block_cell_size);

        /*
         * The number of total value entries per block.
         * Each block holds indices for <m_keys_per_ctrl_block> keys
         * and the values for thes keys may not be more than <max_block_size>
         */
        size_t max_block_size = (1 << sizeof(KeyMap_t)*8); //  This is for block
        size_t values_size= 0;
        ValueEntry* m_map = nullptr;

        // The packed layout of a query run (PackedLayout, Pack): m_packed holds the values, m_map is
        // freed. The key map is the same; a key's values start at bit (its first slot) * m_slot_bits.
        uint8_t* m_packed = nullptr;
        uint64_t m_packed_bytes = 0;
        uint64_t m_packed_entries = 0;
        unsigned m_taxid_bits = 0, m_gene_bits = 0, m_pos_bits = 0, m_entry_bits = 0, m_slot_bits = 0;
        uint64_t m_taxid_mask = 0, m_gene_mask = 0, m_pos_mask = 0, m_entry_mask = 0;

        // What this build writes (plus kFeatureReferenceFingerprint once SetReferenceFingerprint is
        // called); replaced by the index's own features on Load.
        uint32_t m_features = kFeatureFullSyncmerMask | kFeatureCorrectUniqueTwoFlags | kFeatureCheckedSingleEntries;
        ReferenceFingerprint m_reference{};

        uint64_t m_found_counter = 0;

        void Print() {
            std::cout << "m_keys_per_ctrl_block:               " << m_keys_per_ctrl_block << std::endl;
            std::cout << "ctrl_block_frequency_bitshift:       " << ctrl_block_frequency_bitshift << std::endl;
            std::cout << "ctrl_block_byte_size:                " << ctrl_block_byte_size << std::endl;
            std::cout << "ctrl_block_cell_size:                " << ctrl_block_cell_size << std::endl;
            std::cout << "ctrl_block_byte_size_shift:          " << ctrl_block_byte_size_shift << std::endl;
            std::cout << "ctrl_block_key_mask:                 " << ctrl_block_key_mask << std::endl;
            std::cout << "keymap_size_total:                   " << keymap_size_total << std::endl;
            std::cout << "max_block_size:                      " << max_block_size << std::endl;
            std::cout << "values_size:                         " << values_size << std::endl;
        }

        Seedmap(Seedmap const& other) = delete;

        // The key map is counted up in place while building, so it must start zeroed. calloc gives
        // zeroed memory without touching it (large blocks are fresh zero pages), which keeps the
        // default-constructed map of a query run, replaced by Load, cheap.
        void AllocateKeymap(size_t cells) {
            std::free(m_keymap);
            m_keymap = static_cast<KeyMap_t*>(std::calloc(cells, sizeof(KeyMap_t)));
            if (!m_keymap) {
                std::cerr << "Cannot allocate the index key map (" << cells * sizeof(KeyMap_t) << " bytes)" << std::endl;
                exit(8);
            }
            AdviseHugePages(m_keymap, cells * sizeof(KeyMap_t));
        }

        // Every k-mer lookup reads the key map (~3 GB) and the values at random, so with 4 KB pages
        // nearly each one also misses the TLB. Transparent huge pages, which Linux gives to memory
        // that asks for them in its default "madvise" mode, cut seeding time by about a third.
        // Called before the memory is first touched (calloc and malloc hand out untouched pages for
        // blocks this large), so its pages fault in as huge pages. Without THP this does nothing.
        static void AdviseHugePages(void* data, size_t bytes) {
#ifdef MADV_HUGEPAGE
            constexpr size_t kMinBytes = size_t{64} << 20;  // smaller blocks may share heap pages
            if (bytes < kMinBytes) return;
            auto const page = static_cast<uintptr_t>(sysconf(_SC_PAGESIZE));
            auto const begin = (reinterpret_cast<uintptr_t>(data) + page - 1) & ~(page - 1);
            auto const end = (reinterpret_cast<uintptr_t>(data) + bytes) & ~(page - 1);
            if (end > begin) madvise(reinterpret_cast<void*>(begin), end - begin, MADV_HUGEPAGE);
#endif
        }

        // The value array is malloc'd (a loaded index overwrites it anyway; the build zeroes it with
        // calloc), so that loading threads touch its pages first, in parallel. ValueEntry is
        // trivially copyable, so the bytes read into it are its values.
        void AllocateValues(size_t count, bool zero = false) {
            std::free(m_map);
            size_t const bytes = std::max<size_t>(count, 1) * sizeof(ValueEntry);
            m_map = static_cast<ValueEntry*>(zero ? std::calloc(std::max<size_t>(count, 1), sizeof(ValueEntry)) : std::malloc(bytes));
            if (!m_map) {
                std::cerr << "Cannot allocate the index values (" << bytes << " bytes)" << std::endl;
                exit(8);
            }
            AdviseHugePages(m_map, bytes);
        }

        [[noreturn]] static void InvalidIndex(std::string const& name, std::string const& reason) {
            std::cerr << "Invalid index " << name << ": " << reason << std::endl;
            exit(8);
        }

    public:
        // How a query run holds the values (Pack; Load with a layout). The file keeps every entry in a
        // 64-bit slot (taxid, gene and position of 20 bits, 2 flags) and a 32-bit flex cell per entry of a
        // key with several values, in ceil(S/3) of the key's S slots. In memory an entry takes
        // taxid_bits + gene_bits + pos_bits + 2 bits, the widths the reference needs (18 + 8 + 14 + 2 at
        // GTDB r226: 42 of the 64), and the flex cells stay 32 bits. A key's region starts at bit
        // (first slot) * SlotBits(): its e flex cells (if S >= 2), then its e entries at EntryBits() each,
        // all at bit offsets, every entry still addressable on its own. The key map is unchanged (its
        // offsets count the file's slots), so the regions of a load's chunks are known before they are
        // decoded. SlotBits() = ceil(2 (32 + W) / 3) (at least W) holds every S: a key of S slots has
        // e = S - ceil(S/3) <= 2S/3 entries. 50 bits per slot at r226: 27.3 GB for the 35 GB of slots
        // (docs/claude/2026-10-03-memory-audit).
        struct PackedLayout {
            unsigned taxid_bits = 20, gene_bits = 20, pos_bits = 20;

            unsigned EntryBits() const { return taxid_bits + gene_bits + pos_bits + 2; }

            unsigned SlotBits() const {
                unsigned const w = EntryBits();
                return std::max(w, (2 * (32 + w) + 2) / 3);
            }

            // Bits for values up to `max`: 1 to 20 (the file's field width).
            static unsigned Bits(uint64_t max) {
                unsigned b = 1;
                while (b < SEEDMAP_TAXID_BITS && (max >> b) != 0) b++;
                return b;
            }

            // The layout for a reference whose taxids, gene ids and positions go up to these values.
            static PackedLayout For(uint64_t max_taxid, uint64_t max_gene, uint64_t max_position) {
                return { Bits(max_taxid), Bits(max_gene), Bits(max_position) };
            }

            bool operator==(PackedLayout const&) const = default;
        };

        // A key's values in the packed layout (GetPacked).
        struct PackedBlock {
            uint8_t const* entries = nullptr;  // entry i: EntryBits() bits at bit entry_shift + i * EntryBits() from here
            uint8_t const* flex = nullptr;     // flex cell i: 32 bits at bit flex_shift + 32 i from here; nullptr if the key has one value
            uint32_t entry_shift = 0, flex_shift = 0;  // 0-7
            uint32_t size = 0;                 // entries
        };

        Seedmap() {
            AllocateKeymap(keymap_size_total);
        }

        bool IsPacked() const {
            return m_packed != nullptr;
        }

        PackedLayout Layout() const {
            return { m_taxid_bits, m_gene_bits, m_pos_bits };
        }

        uint64_t PackedEntries() const {
            return m_packed_entries;
        }

        // The memory the index takes: the values (packed: entries and flex cells; else the slots) and the key map.
        std::string MemoryDescription() const {
            auto gb = [](uint64_t bytes) {
                char buf[32];
                snprintf(buf, sizeof(buf), "%.1f GB", bytes / 1e9);
                return std::string(buf);
            };
            std::string s;
            if (IsPacked()) {
                s = std::to_string(m_packed_entries) + " entries of " + std::to_string(m_entry_bits) + " bits (taxid " +
                    std::to_string(m_taxid_bits) + ", gene " + std::to_string(m_gene_bits) + ", position " + std::to_string(m_pos_bits) +
                    ", 2 flags) and their 32-bit flex cells, " + std::to_string(m_slot_bits) + " bits per slot of the file: " + gb(m_packed_bytes);
            } else {
                s = std::to_string(values_size) + " slots of 8 bytes: " + gb(values_size * sizeof(ValueEntry));
            }
            return s + "; key map " + gb(keymap_size_total * sizeof(KeyMap_t));
        }

        // Flex cell i of a key's values.
        static inline uint32_t FlexCell(PackedBlock const& block, uint32_t i) {
            uint64_t x;
            std::memcpy(&x, block.flex + 4 * static_cast<size_t>(i), 8);
            return static_cast<uint32_t>(x >> block.flex_shift);
        }

        // Entry i of a key's values, as the bits of a ValueEntry (taxid << 40 | gene << 20 | position, flags at 60).
        inline uint64_t EntryValue(PackedBlock const& block, uint32_t i) const {
            uint64_t const bit = block.entry_shift + static_cast<uint64_t>(i) * m_entry_bits;
            __uint128_t x;
            std::memcpy(&x, block.entries + (bit >> 3), 16);
            uint64_t const v = static_cast<uint64_t>(x >> (bit & 7)) & m_entry_mask;
            uint64_t const taxid = (v >> (m_pos_bits + m_gene_bits)) & m_taxid_mask;
            uint64_t const gene = (v >> m_pos_bits) & m_gene_mask;
            uint64_t const flags = v >> (m_pos_bits + m_gene_bits + m_taxid_bits);
            return taxid << 40 | gene << 20 | (v & m_pos_mask) | flags << 60;
        }

        // The slots of a key's values in the file's layout: the first and how many; false if none.
        inline bool Locate(uint64_t main_key, uint64_t& start, uint64_t& slots) const {
            uint64_t const block_start_idx = ControlBlockIndex(main_key);
            uint64_t const block_end_idx = block_start_idx + m_keys_per_ctrl_block + ctrl_block_cell_size;
            uint64_t block_value_start_idx, block_value_end_idx;
            std::memcpy(&block_value_start_idx, m_keymap + block_start_idx, 8);
            std::memcpy(&block_value_end_idx, m_keymap + block_end_idx, 8);
            if (block_value_end_idx == block_value_start_idx) return false;
            size_t const key_index = block_start_idx + ctrl_block_cell_size + BlockKey(main_key);
            start = block_value_start_idx + m_keymap[key_index];
            uint64_t const end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx + m_keymap[key_index + 1];
            if (end <= start) return false;
            slots = end - start;
            return true;
        }

        // The values of a k-mer's core in the packed layout; false if it has none.
        inline bool GetPacked(uint64_t key, PackedBlock& block) const {
            uint64_t start, slots;
            if (!Locate(MainKey(key), start, slots)) return false;
            BlockOfSlots(start, slots, block);
            return true;
        }

        // The packed values of the key whose slots (in the file's layout) are [start, start + slots), slots >= 1.
        inline void BlockOfSlots(uint64_t start, uint64_t slots, PackedBlock& block) const {
            uint64_t const bit = start * m_slot_bits;
            if (slots >= m_flex_threshold) {
                uint64_t const entries = slots - FlexBlockSize(slots);
                block.size = static_cast<uint32_t>(entries);
                block.flex = m_packed + (bit >> 3);
                block.flex_shift = static_cast<uint32_t>(bit & 7);
                uint64_t const entry_bit = bit + 32 * entries;
                block.entries = m_packed + (entry_bit >> 3);
                block.entry_shift = static_cast<uint32_t>(entry_bit & 7);
            } else {
                block.size = static_cast<uint32_t>(slots);
                block.flex = nullptr;
                block.entries = m_packed + (bit >> 3);
                block.entry_shift = static_cast<uint32_t>(bit & 7);
            }
        }

        // The bit of flag `flag` (0: unique, 1: unique at distance two) of entry i of a key's packed values.
        uint64_t FlagBit(PackedBlock const& block, uint32_t i, unsigned flag) const {
            return static_cast<uint64_t>(block.entries - m_packed) * 8 + block.entry_shift + static_cast<uint64_t>(i) * m_entry_bits +
                   m_taxid_bits + m_gene_bits + m_pos_bits + flag;
        }

        // The flags of entry i of a key's packed values, changed atomically a byte at a time: --build's
        // uniqueness check clears them from any thread, the unique k-mer statistics set them by block ranges,
        // and an entry's bytes may be shared with a neighbouring key's.
        void ClearUniqueFlag(PackedBlock const& block, uint32_t i) {
            uint64_t const bit = FlagBit(block, i, 0);
            std::atomic_ref<uint8_t>(m_packed[bit >> 3]).fetch_and(static_cast<uint8_t>(~(1u << (bit & 7))), std::memory_order_relaxed);
        }

        void SetUniqueDistanceTwoFlag(PackedBlock const& block, uint32_t i) {
            uint64_t const bit = FlagBit(block, i, 1);
            std::atomic_ref<uint8_t>(m_packed[bit >> 3]).fetch_or(static_cast<uint8_t>(1u << (bit & 7)), std::memory_order_relaxed);
        }

        // The cells [first, first + n) of the values in the file's 8-byte layout, from the packed values: a
        // key's flex cells two to a slot (the second 0 for an odd number of values), then its entries.
        void UnpackSlots(uint64_t first, uint64_t n, uint64_t* out) const {
            if (n == 0) return;
            uint64_t const cpb = m_keys_per_ctrl_block + ctrl_block_cell_size, blocks = keymap_size / m_keys_per_ctrl_block;
            // The block holding slot `first`: the last one starting at or before it (empty blocks start where the next does).
            uint64_t lo = 0, hi = blocks;
            while (hi - lo > 1) {
                uint64_t const mid = lo + (hi - lo) / 2;
                (BlockStart(mid) <= first ? lo : hi) = mid;
            }
            uint64_t const last = first + n;
            for (uint64_t block = lo; block < blocks && BlockStart(block) < last; block++) {
                uint64_t const begin = BlockStart(block), end = BlockStart(block + 1);
                if (end <= first) continue;
                for (size_t j = 0; j < m_keys_per_ctrl_block; j++) {
                    uint64_t const start = begin + m_keymap[block * cpb + ctrl_block_cell_size + j];
                    uint64_t const stop = j + 1 == m_keys_per_ctrl_block ? end : begin + m_keymap[block * cpb + ctrl_block_cell_size + j + 1];
                    if (stop <= start || stop <= first || start >= last) continue;
                    uint64_t const S = stop - start;
                    PackedBlock values;
                    BlockOfSlots(start, S, values);
                    uint64_t const f = S >= m_flex_threshold ? FlexBlockSize(S) : 0;
                    for (uint64_t s = std::max(start, first); s < std::min(stop, last); s++) out[s - first] = SlotCell(values, f, s - start);
                }
            }
        }

        // The cells of the key whose slots are [start, start + S) (S >= 1), as UnpackSlots gives them,
        // without looking for the key.
        void UnpackKey(uint64_t start, uint64_t S, uint64_t* out) const {
            PackedBlock values;
            BlockOfSlots(start, S, values);
            uint64_t const f = S >= m_flex_threshold ? FlexBlockSize(S) : 0;
            for (uint64_t j = 0; j < S; j++) out[j] = SlotCell(values, f, j);
        }

        // Slot j of a key's values (f of its slots are flex slots) in the file's 8-byte layout.
        uint64_t SlotCell(PackedBlock const& values, uint64_t f, uint64_t j) const {
            if (j >= f) return EntryValue(values, static_cast<uint32_t>(j - f));
            uint64_t const a = 2 * j, b = a + 1;
            return FlexCell(values, static_cast<uint32_t>(a)) |
                   (b < values.size ? uint64_t{FlexCell(values, static_cast<uint32_t>(b))} << 32 : 0);
        }

        // ORs the low w bits of x into the zeroed bit array `base` at `bit`. Bytes in [safe_begin, safe_end)
        // are this thread's alone (a plain read-modify-write of 8 bytes); any other byte may be shared
        // with the range a neighbouring thread packs, so it is OR'd a byte at a time, atomically.
        static void PutBits(uint8_t* base, uint64_t bit, uint64_t x, unsigned w, uint64_t safe_begin, uint64_t safe_end) {
            uint64_t const byte = bit >> 3;
            unsigned const shift = static_cast<unsigned>(bit & 7);
            if (shift + w <= 64 && byte >= safe_begin && byte + 8 <= safe_end) {
                uint64_t cur;
                std::memcpy(&cur, base + byte, 8);
                cur |= x << shift;
                std::memcpy(base + byte, &cur, 8);
                return;
            }
            __uint128_t const v = static_cast<__uint128_t>(x) << shift;
            for (unsigned k = 0; k < (shift + w + 7) / 8; k++) {
                auto const b = static_cast<uint8_t>(v >> (8 * k));
                if (b) std::atomic_ref<uint8_t>(base[byte + k]).fetch_or(b, std::memory_order_relaxed);
            }
        }

        // The w bits (w <= 64) at `bit` of the bit array `base`, as PutBits writes them: bytes outside
        // [safe_begin, safe_end) are read atomically, as another thread may be OR'ing into them.
        static uint64_t GetBits(uint8_t* base, uint64_t bit, unsigned w, uint64_t safe_begin, uint64_t safe_end) {
            uint64_t const byte = bit >> 3;
            unsigned const shift = static_cast<unsigned>(bit & 7);
            unsigned const bytes = (shift + w + 7) / 8;
            __uint128_t x = 0;
            if (byte >= safe_begin && byte + bytes <= safe_end) {
                std::memcpy(&x, base + byte, bytes);
            } else {
                for (unsigned k = 0; k < bytes; k++) {
                    x |= static_cast<__uint128_t>(std::atomic_ref<uint8_t>(base[byte + k]).load(std::memory_order_relaxed)) << (8 * k);
                }
            }
            uint64_t const mask = w == 64 ? ~uint64_t{0} : (uint64_t{1} << w) - 1;
            return static_cast<uint64_t>(x >> shift) & mask;
        }

        // A value (the bits of a ValueEntry: taxid << 40 | gene << 20 | position, flags at 60) in the packed
        // layout's fields; false if a field does not fit them.
        bool PackValue(uint64_t v, uint64_t& packed) const {
            uint64_t const taxid = (v >> 40) & 0xfffff, gene = (v >> 20) & 0xfffff, pos = v & 0xfffff, flags = (v >> 60) & 3;
            if (taxid > m_taxid_mask || gene > m_gene_mask || pos > m_pos_mask) return false;
            packed = taxid << (m_gene_bits + m_pos_bits) | gene << m_pos_bits | pos | flags << (m_taxid_bits + m_gene_bits + m_pos_bits);
            return true;
        }

        // Converts the loaded 8-byte layout into the packed one with `threads` threads and frees it.
        // Exits if a value's fields do not fit the layout (the index was built against another reference).
        void Pack(PackedLayout const& layout, int threads) {
            SetLayout(layout);
            AllocatePacked();
            uint64_t const blocks = keymap_size / m_keys_per_ctrl_block;
            size_t const parts = static_cast<size_t>(std::max(threads, 1)) * 16;
            std::atomic<uint64_t> entries{0};
            std::string const error = zstd::ParallelFor(parts, threads, [&](size_t part, size_t) -> std::string {
                uint64_t const b0 = blocks * part / parts, b1 = blocks * (part + 1) / parts;
                if (b0 == b1) return "";
                uint64_t bad = 0;
                auto const n = PackBlocks(b0, b1, reinterpret_cast<uint64_t const*>(m_map), 0, BlockStart(b0), BlockStart(b1), bad);
                if (!n) return BadValueMessage(bad);
                entries += *n;
                return "";
            });
            if (!error.empty()) InvalidIndex("index.prx", error);
            m_packed_entries = entries;
            std::free(m_map);
            m_map = nullptr;
        }

    private:
        void SetLayout(PackedLayout const& layout) {
            m_taxid_bits = layout.taxid_bits;
            m_gene_bits = layout.gene_bits;
            m_pos_bits = layout.pos_bits;
            m_entry_bits = layout.EntryBits();
            m_slot_bits = layout.SlotBits();
            m_taxid_mask = (uint64_t{1} << m_taxid_bits) - 1;
            m_gene_mask = (uint64_t{1} << m_gene_bits) - 1;
            m_pos_mask = (uint64_t{1} << m_pos_bits) - 1;
            m_entry_mask = (uint64_t{1} << m_entry_bits) - 1;
        }

        // Zeroed (the packers OR their bits in) and untouched until a thread writes its range, as the
        // genes' arena (GenomeLoader::LoadAllGenomes); 16 bytes more, read past the last entry.
        void AllocatePacked() {
            std::free(m_packed);
            m_packed_bytes = (values_size * m_slot_bits + 7) / 8;
            m_packed = static_cast<uint8_t*>(std::calloc(m_packed_bytes + 16, 1));
            if (!m_packed) {
                std::cerr << "Cannot allocate the index values (" << m_packed_bytes + 16 << " bytes)" << std::endl;
                exit(8);
            }
            AdviseHugePages(m_packed, m_packed_bytes + 16);
        }

        // The first slot of control block `block` (block == number of blocks: the number of slots).
        uint64_t BlockStart(uint64_t block) const {
            uint64_t v;
            std::memcpy(&v, m_keymap + block * (m_keys_per_ctrl_block + ctrl_block_cell_size), 8);
            return v;
        }

        static std::string BadValueMessage(uint64_t value) {
            ValueEntry entry;
            entry.value = value;
            auto const [taxid, gene, pos] = entry.Get();
            return "holds a value outside the reference's ranges (taxid " + std::to_string(taxid) + ", gene " +
                   std::to_string(gene) + ", position " + std::to_string(pos) + "): it was built against a different reference";
        }

        // Packs one key's S slots (at `cells`, the file's layout) into its region at bit slot * m_slot_bits.
        bool PackKey(uint64_t slot, uint64_t S, uint64_t const* cells, uint64_t safe_begin, uint64_t safe_end, uint64_t& bad) {
            uint64_t bit = slot * m_slot_bits;
            uint64_t e = S;
            uint64_t const* entries = cells;
            if (S >= m_flex_threshold) {
                uint64_t const f = FlexBlockSize(S);
                e = S - f;
                entries = cells + f;
                auto const* flex = reinterpret_cast<uint32_t const*>(cells);
                for (uint64_t i = 0; i < e; i++) PutBits(m_packed, bit + 32 * i, flex[i], 32, safe_begin, safe_end);
                bit += 32 * e;
            }
            for (uint64_t i = 0; i < e; i++, bit += m_entry_bits) {
                uint64_t packed = 0;
                if (!PackValue(entries[i], packed)) {
                    bad = entries[i];
                    return false;
                }
                PutBits(m_packed, bit, packed, m_entry_bits, safe_begin, safe_end);
            }
            return true;
        }

        // Packs the keys of control blocks [first_block, end_block), whose stored values lie at
        // values[slot - values_first_slot] and whose slots are [range_first_slot, range_end_slot) (the
        // end also bounds the last block, so nothing of the next range is read). The bytes of the
        // region shared with the neighbouring ranges are OR'd atomically (PutBits). Returns the entries
        // packed, or nullopt with `bad` the first value whose fields do not fit the layout.
        std::optional<uint64_t> PackBlocks(uint64_t first_block, uint64_t end_block, uint64_t const* values, uint64_t values_first_slot,
                                           uint64_t range_first_slot, uint64_t range_end_slot, uint64_t& bad) {
            uint64_t const first_bit = range_first_slot * m_slot_bits, end_bit = range_end_slot * m_slot_bits;
            uint64_t const safe_begin = (first_bit >> 3) + ((first_bit & 7) ? 1 : 0);
            uint64_t const safe_end = end_bit >> 3;
            uint64_t const cpb = m_keys_per_ctrl_block + ctrl_block_cell_size;
            uint64_t entries = 0;
            for (uint64_t block = first_block; block < end_block; block++) {
                uint64_t const ctrl = block * cpb;
                uint64_t const begin = BlockStart(block);
                uint64_t const end = block + 1 == end_block ? range_end_slot : BlockStart(block + 1);
                if (begin == end) continue;
                for (size_t j = 0; j < m_keys_per_ctrl_block; j++) {
                    uint64_t const start = begin + m_keymap[ctrl + ctrl_block_cell_size + j];
                    uint64_t const stop = j + 1 == m_keys_per_ctrl_block ? end : begin + m_keymap[ctrl + ctrl_block_cell_size + j + 1];
                    if (stop <= start) continue;
                    uint64_t const S = stop - start;
                    if (!PackKey(start, S, values + (start - values_first_slot), safe_begin, safe_end, bad)) return std::nullopt;
                    entries += S >= m_flex_threshold ? S - FlexBlockSize(S) : S;
                }
            }
            return entries;
        }

    public:

        Seedmap(std::string file) {
            Load(file);
        };

        ~Seedmap() {
            std::free(m_keymap);
            std::free(m_map);
            std::free(m_packed);
        }

        inline uint64_t KeymapIndex(uint64_t key) const {
            return key + ((key >> ctrl_block_frequency_bitshift) << ctrl_block_byte_size_shift) + ctrl_block_cell_size;
        }

        inline uint64_t ControlBlockIndex(uint64_t key) const {
            return ((key >> ctrl_block_frequency_bitshift) << ctrl_block_frequency_bitshift) + ((key >> ctrl_block_frequency_bitshift) << (ctrl_block_byte_size_shift));
        }

        inline uint64_t BlockKey(uint64_t key) const {
            return key & ctrl_block_key_mask;
        }

        inline static uint64_t DivisionByTwoCeiling(uint64_t x) {
            return (x & 1llu) ? (x >> 1) + 1 : (x >> 1);
        }

        bool UsesFullSyncmerMask() const {
            return m_features & kFeatureFullSyncmerMask;
        }

        std::string FeatureDescription() const {
            std::string description = (m_features & kFeatureFullSyncmerMask) ? "full s-mer syncmers" : "legacy 4-base syncmers";
            description += (m_features & kFeatureCorrectUniqueTwoFlags) ? ", per-entry unique-distance flags"
                                                                        : ", legacy unique-distance flags";
            description += HasReferenceFingerprint() ? ", reference fingerprint" : ", no reference fingerprint";
            description += ChecksSingleEntryUniques() ? ", checked single-entry uniques"
                                                      : ", unchecked single-entry uniques (rebuild for correct unique_kmers.tsv)";
            return description;
        }

        // Whether k-mer cores that occur once in the index were checked against the full reference
        // (before, they all stayed unique, which made every gene hittable).
        bool ChecksSingleEntryUniques() const {
            return m_features & kFeatureCheckedSingleEntries;
        }

        void SetReferenceFingerprint(ReferenceFingerprint const& fingerprint) {
            m_reference = fingerprint;
            m_features |= kFeatureReferenceFingerprint;
        }

        bool HasReferenceFingerprint() const {
            return m_features & kFeatureReferenceFingerprint;
        }

        ReferenceFingerprint const& GetReferenceFingerprint() const {
            return m_reference;
        }

        void SaveHeader(std::ostream& ofs) {
            const uint64_t file_magic = kFileMagicV2;
            const uint32_t version_major = kIndexVersionMajor;
            const uint32_t version_minor = kIndexVersionMinor;
            const uint32_t version_patch = kIndexVersionPatch;
            const uint32_t features = m_features;
            ofs.write((char *) &file_magic, sizeof(file_magic));
            ofs.write((char *) &version_major, sizeof(version_major));
            ofs.write((char *) &version_minor, sizeof(version_minor));
            ofs.write((char *) &version_patch, sizeof(version_patch));
            ofs.write((char *) &features, sizeof(features));
            if (HasReferenceFingerprint()) {
                ofs.write((char *) &m_reference.map_hash, sizeof(m_reference.map_hash));
                ofs.write((char *) &m_reference.fna_size, sizeof(m_reference.fna_size));
            }
        }

        void LoadHeader(std::istream& ifs) {
            uint64_t file_magic = 0;
            uint32_t version_major = 0;
            uint32_t version_minor = 0;
            uint32_t version_patch = 0;

            ifs.read((char *) &file_magic, sizeof(file_magic));
            ifs.read((char *) &version_major, sizeof(version_major));
            ifs.read((char *) &version_minor, sizeof(version_minor));
            ifs.read((char *) &version_patch, sizeof(version_patch));

            if (!ifs) {
                std::cerr << "Failed to read index header from index.prx" << std::endl;
                exit(8);
            }

            m_features = 0;
            if (file_magic == kFileMagicV2) {
                ifs.read((char *) &m_features, sizeof(m_features));
                if (!ifs) {
                    std::cerr << "Failed to read index header from index.prx" << std::endl;
                    exit(8);
                }
                if (m_features & ~kKnownFeatures) {
                    std::cerr << "index.prx uses features this protal does not know (" << m_features
                              << "); it was built by a newer protal" << std::endl;
                    exit(8);
                }
                if (HasReferenceFingerprint()) {
                    ifs.read((char *) &m_reference.map_hash, sizeof(m_reference.map_hash));
                    ifs.read((char *) &m_reference.fna_size, sizeof(m_reference.fna_size));
                    if (!ifs) {
                        std::cerr << "Failed to read index header from index.prx" << std::endl;
                        exit(8);
                    }
                }
            } else if (file_magic != kFileMagic) {
                std::cerr << "Unsupported index.prx format: missing or invalid file header. "
                          << "This binary expects indices written by protal version "
                          << kOldestCompatibleIndexVersionMajor << "."
                          << kOldestCompatibleIndexVersionMinor << "."
                          << kOldestCompatibleIndexVersionPatch
                          << " or newer" << std::endl;
                exit(8);
            }

            const bool too_old =
                    version_major < kOldestCompatibleIndexVersionMajor ||
                    (version_major == kOldestCompatibleIndexVersionMajor &&
                     version_minor < kOldestCompatibleIndexVersionMinor) ||
                    (version_major == kOldestCompatibleIndexVersionMajor &&
                     version_minor == kOldestCompatibleIndexVersionMinor &&
                     version_patch < kOldestCompatibleIndexVersionPatch);

            const bool too_new =
                    version_major > kIndexVersionMajor ||
                    (version_major == kIndexVersionMajor &&
                     version_minor > kIndexVersionMinor) ||
                    (version_major == kIndexVersionMajor &&
                     version_minor == kIndexVersionMinor &&
                     version_patch > kIndexVersionPatch);

            if (too_old || too_new) {
                std::cerr << "Unsupported index.prx version "
                          << version_major << "." << version_minor << "." << version_patch
                          << ". This binary supports indices written by protal versions "
                          << kOldestCompatibleIndexVersionMajor << "."
                          << kOldestCompatibleIndexVersionMinor << "."
                          << kOldestCompatibleIndexVersionPatch
                          << " to "
                          << kIndexVersionMajor << "."
                          << kIndexVersionMinor << "."
                          << kIndexVersionPatch << std::endl;
                exit(8);
            }
        }

        void Save(std::string file) {
            std::ofstream ofs(file, std::ios::binary);
            Save(ofs);
            ofs.close();
            if (ofs.fail()) {
                std::cerr << "Writing the index " << file << " failed" << std::endl;
                exit(8);
            }
        }

        // Bytes Save writes: file header, six layout fields, key map, values.
        uint64_t SerializedSize() const {
            uint64_t const fingerprint = HasReferenceFingerprint() ? sizeof(m_reference.map_hash) + sizeof(m_reference.fna_size) : 0;
            return sizeof(uint64_t) + 4 * sizeof(uint32_t) + fingerprint + 6 * sizeof(size_t) +
                   keymap_size_total * sizeof(KeyMap_t) + values_size * sizeof(ValueEntry);
        }

        // The bytes Save writes before the key map: file header and the six layout fields.
        std::string HeaderBytes() {
            std::ostringstream ofs(std::ios::binary);
            SaveHeader(ofs);
            ofs.write((char *) &keymap_size, sizeof(keymap_size));
            ofs.write((char *) &keymap_size_total, sizeof(keymap_size_total));
            ofs.write((char *) &ctrl_block_byte_size, sizeof(ctrl_block_byte_size));
            ofs.write((char *) &m_keys_per_ctrl_block, sizeof(m_keys_per_ctrl_block));
            ofs.write((char *) &ctrl_block_frequency_bitshift, sizeof(ctrl_block_frequency_bitshift));
            ofs.write((char *) &values_size, sizeof(values_size));
            return ofs.str();
        }

        // The raw index (header, key map, values in the 8-byte layout); packed values are written unpacked, in batches.
        void Save(std::ostream& ofs) {
//            SortForKeys();
            std::string const header = HeaderBytes();
            ofs.write(header.data(), static_cast<std::streamsize>(header.size()));
            ofs.write((char *) m_keymap, sizeof(*m_keymap) * keymap_size_total);
            if (!IsPacked()) {
                ofs.write((char *) m_map, sizeof(*m_map) * values_size);
                return;
            }
            constexpr uint64_t kBatch = uint64_t{1} << 20;
            std::vector<uint64_t> cells(kBatch);
            for (uint64_t first = 0; first < values_size && ofs; first += kBatch) {
                uint64_t const n = std::min(kBatch, values_size - first);
                UnpackSlots(first, n, cells.data());
                ofs.write(reinterpret_cast<char const*>(cells.data()), static_cast<std::streamsize>(n * sizeof(uint64_t)));
            }
        }

        index_codec::Layout CodecLayout() const {
            return { keymap_size / m_keys_per_ctrl_block, m_keys_per_ctrl_block, ctrl_block_cell_size, values_size };
        }

        // The values as the column format's writer and Verify read them: the 8-byte layout in place, or the packed
        // values unpacked a key at a time (UnpackKey; UnpackSlots for other ranges).
        index_codec::Cells ValueCells() const {
            if (!IsPacked()) return { reinterpret_cast<uint64_t const*>(m_map) };
            return { this,
                     [](void const* map, uint64_t first, uint64_t n, uint64_t* out) { static_cast<Seedmap const*>(map)->UnpackKey(first, n, out); },
                     [](void const* map, uint64_t first, uint64_t n, uint64_t* out) { static_cast<Seedmap const*>(map)->UnpackSlots(first, n, out); } };
        }

        // Writes the index compressed in columns (IndexCodec.h), then reads the file back chunk by
        // chunk and compares it with the index. Returns an error message, empty on success.
        std::string SaveCompressed(std::string const& path, zstd::Params const& params, uint64_t& written, size_t& raw_chunks) {
            std::string error;
            std::string const header = HeaderBytes();
            index_codec::Cells const values = ValueCells();
            auto const bytes = index_codec::Write(path, header, CodecLayout(), m_keymap, values, params, error, &raw_chunks);
            if (!bytes) return error;
            written = *bytes;
            error = index_codec::Verify(path, header, CodecLayout(), m_keymap, values, params.threads);
            return error.empty() ? error : "reading it back: " + error;
        }

        // Reads an index written by Save. Every size in the file is checked against the layout this
        // build uses and against the bytes the stream actually holds, so a truncated, corrupt or
        // incompatible index stops here with a message instead of causing out-of-bounds reads later.
        void Load(std::istream &ifs, std::string const& name = "index.prx") {
            uint64_t const expected_bytes = LoadLayout(ifs, name);

            // The rest of the file must hold exactly the key map and the values (when the size is known).
            auto const data_start = ifs.tellg();
            if (data_start != std::streampos(-1) && ifs.seekg(0, std::ios::end)) {
                auto const file_end = ifs.tellg();
                ifs.seekg(data_start);
                if (file_end != std::streampos(-1) && static_cast<size_t>(file_end - data_start) != expected_bytes) {
                    InvalidIndex(name, std::to_string(static_cast<size_t>(file_end - data_start)) + " bytes of data, expected " +
                                       std::to_string(expected_bytes) + " (truncated or corrupt file?)");
                }
            }
            ifs.clear();

            AllocateKeymap(keymap_size_total);
            ifs.read((char *) m_keymap, sizeof(*m_keymap) * (keymap_size_total));
            AllocateValues(values_size);
            ifs.read((char *) m_map, sizeof(*m_map) * (values_size));
            if (ifs.bad()) InvalidIndex(name, "the file could not be read or decompressed (truncated or corrupt file?)");
            if (!ifs) InvalidIndex(name, "the file ends before all index data was read (truncated or corrupt file?)");
            // Streams whose size is unknown (zstd) are checked here. Reaching the end also makes
            // zstd verify the frame's content checksum.
            if (!std::char_traits<char>::eq_int_type(ifs.peek(), std::char_traits<char>::eof())) {
                InvalidIndex(name, "unexpected data after the index (corrupt file?)");
            }
            if (ifs.bad()) InvalidIndex(name, "the file could not be read or decompressed (truncated or corrupt file?)");
        }

        // Reads index.prx or index.prx.zst with `threads` threads: the column format --build writes
        // (IndexCodec.h) chunk by chunk, a raw file or a seekable zstd file of the raw bytes straight
        // into the key map and values, any other zstd file in one stream. With a layout the values
        // end up packed (PackedLayout): the column format is packed as it is decoded, the other
        // formats after they are read (their 8-byte layout is in memory until then).
        void Load(std::string file, int threads = 1, PackedLayout const* pack = nullptr) {
            zstd::InputFile in(file);
            if (!in.IsOpen()) InvalidIndex(file, "cannot open the file");
            if (in.Compressed()) {
                std::string error;
                auto const table = zstd::ReadSeekTable(file, error);
                if (!error.empty()) InvalidIndex(file, error + " (truncated or corrupt file?)");
                if (!table) {
                    if (index_codec::StartsWithMagic(file)) {
                        InvalidIndex(file, "the seek table at its end is missing (truncated or corrupt file?)");
                    }
                    Load(in.Stream(), file);
                    if (pack) Pack(*pack, threads);
                    return;
                }
                auto const container = index_codec::ReadContainer(file, *table, error);
                if (!error.empty()) InvalidIndex(file, error + " (truncated or corrupt file?)");
                if (container) {
                    LoadColumns(file, file, *table, *container, threads, pack);
                    return;
                }
            }
            uint64_t const data_bytes = LoadLayout(in.Stream(), file);
            auto const position = in.Stream().tellg();
            if (position == std::streampos(-1)) InvalidIndex(file, "cannot read the file");
            uint64_t const header_bytes = static_cast<uint64_t>(position);
            // Check the size first (raw: file size, seekable: its seek table), before allocating.
            if (auto const total = zstd::UncompressedSize(file); total && *total != header_bytes + data_bytes) {
                uint64_t const found = *total > header_bytes ? *total - header_bytes : 0;
                InvalidIndex(file, std::to_string(found) + " bytes of data, expected " + std::to_string(data_bytes) +
                                   " (truncated or corrupt file?)");
            }
            AllocateKeymap(keymap_size_total);
            AllocateValues(values_size);
            uint64_t const keymap_bytes = sizeof(*m_keymap) * keymap_size_total;
            IndexSink sink(header_bytes, reinterpret_cast<char*>(m_keymap), keymap_bytes, reinterpret_cast<char*>(m_map),
                           sizeof(*m_map) * values_size);
            std::string error;
            uint64_t const delivered = zstd::ParallelRead(file, threads, sink, error);
            if (!error.empty()) InvalidIndex(file, error + " (truncated or corrupt file?)");
            if (delivered != header_bytes + data_bytes) {
                uint64_t const found = delivered > header_bytes ? delivered - header_bytes : 0;
                InvalidIndex(file, std::to_string(found) + " bytes of data, expected " + std::to_string(data_bytes) +
                                   " (truncated or corrupt file?)");
            }
            if (pack) Pack(*pack, threads);
        }

        // The index as a database file: a file on disk (Load(file, threads) above), or the member of a
        // single-file database (Database.h), which holds the index's frames in the column format.
        void Load(db::DbFile const& file, int threads = 1, PackedLayout const* pack = nullptr) {
            if (!file.InBundle()) {
                Load(file.Path(), threads, pack);
                return;
            }
            if (!file.Exists()) InvalidIndex(file.Name(), "the database has no index");
            std::string error;
            auto const container = index_codec::ReadContainer(file.Path(), file.Frames(), error);
            if (!error.empty()) InvalidIndex(file.Name(), error + " (truncated or corrupt file?)");
            if (!container) InvalidIndex(file.Name(), "not an index in protal's column format (corrupt file?)");
            LoadColumns(file.Path(), file.Name(), file.Frames(), *container, threads, pack);
        }

        // The column format (IndexCodec.h) in the frames `table` lists in the file at path: the raw
        // header from the container, then the chunks, decoded in parallel into the key map and
        // values. name is used in messages. With a layout each chunk's values are packed key by key as
        // they are decoded (PackKey), so the 8-byte layout is never held, not even a chunk of it: the
        // chunk's region in the packed values follows from its first slot.
        void LoadColumns(std::string const& path, std::string const& name, zstd::SeekTable const& table,
                         index_codec::Container const& container, int threads, PackedLayout const* pack = nullptr) {
            std::istringstream header(container.index_header);
            LoadLayout(header, name);
            auto const expected = CodecLayout();
            auto const& l = container.layout;
            if (l.blocks != expected.blocks || l.keys_per_block != expected.keys_per_block ||
                l.ctrl_cells != expected.ctrl_cells || l.values != expected.values) {
                InvalidIndex(name, "the layout of the compressed index does not match its header (corrupt file?)");
            }
            AllocateKeymap(keymap_size_total);
            if (!pack) {
                AllocateValues(values_size);
                std::string const error = index_codec::Decode(path, table, container, m_keymap, reinterpret_cast<uint64_t*>(m_map), threads);
                if (!error.empty()) InvalidIndex(name, error + " (truncated or corrupt file?)");
                return;
            }
            SetLayout(*pack);
            AllocatePacked();
            // Each key's cells are packed as the chunk is decoded (index_codec::ValueSink): no thread holds its chunk's
            // values in the 8-byte layout, ~50 MB each at r226. A chunk stored raw is packed from its cells at once.
            struct Packing {
                Seedmap* map;
                index_codec::Chunk const* chunk;
                uint64_t safe_begin, safe_end;  // the bytes of the chunk's region no other chunk shares (PackBlocks)
                uint64_t entries = 0;
                uint64_t bad = 0;
            };
            index_codec::detail::ValueSink sink;
            sink.key = [](void* context, uint64_t slot, uint64_t const* cells, uint64_t n) {
                auto& p = *static_cast<Packing*>(context);
                if (!p.map->PackKey(slot, n, cells, p.safe_begin, p.safe_end, p.bad)) return false;
                p.entries += n >= p.map->m_flex_threshold ? n - p.map->FlexBlockSize(n) : n;
                return true;
            };
            sink.range = [](void* context, uint64_t first_slot, uint64_t const* cells, uint64_t) {
                auto& p = *static_cast<Packing*>(context);
                auto const& ch = *p.chunk;
                auto const n = p.map->PackBlocks(ch.first_block, ch.first_block + ch.blocks, cells, first_slot, ch.first_value,
                                                 ch.first_value + ch.values, p.bad);
                if (!n) return false;
                p.entries += *n;
                return true;
            };
            std::atomic<uint64_t> entries{0};
            std::string const error = zstd::ForEachFrame(path, table, 1, threads,
                    [&](size_t frame, char const* data, size_t size, size_t) -> std::string {
                index_codec::Chunk const& ch = container.chunks[frame - 1];
                uint64_t const first_bit = ch.first_value * m_slot_bits, end_bit = (ch.first_value + ch.values) * m_slot_bits;
                Packing packing{ this, &ch, (first_bit >> 3) + ((first_bit & 7) ? 1 : 0), end_bit >> 3 };
                auto own = sink;
                own.context = &packing;
                std::string const e = index_codec::detail::DecodeChunk(data, size, l, ch, m_keymap + ch.first_block * l.CellsPerBlock(),
                                                                       nullptr, nullptr, &own);
                if (e == index_codec::detail::kSinkRefused) return BadValueMessage(packing.bad);
                if (!e.empty()) return "chunk " + std::to_string(frame) + " of " + std::to_string(container.chunks.size()) + ": " + e;
                entries += packing.entries;
                return "";
            });
            // A value outside the layout is the index's content (another reference); anything else is the file.
            if (!error.empty()) InvalidIndex(name, error.rfind("holds a value", 0) == 0 ? error : error + " (truncated or corrupt file?)");
            std::memcpy(m_keymap + l.blocks * l.CellsPerBlock(), &l.values, 8);  // final control block
            m_packed_entries = entries;
        }

        // Header and layout fields of an index written by Save; every size is checked against the
        // layout this build uses. Sets up that layout and returns the number of data bytes (key map
        // and values) that must follow.
        uint64_t LoadLayout(std::istream &ifs, std::string const& name) {
            LoadHeader(ifs);
            size_t file_keymap_size = 0, file_keymap_size_total = 0, file_ctrl_block_byte_size = 0;
            size_t file_keys_per_ctrl_block = 0, file_bitshift = 0, file_values_size = 0;
            ifs.read((char *) &file_keymap_size, sizeof(file_keymap_size));
            ifs.read((char *) &file_keymap_size_total, sizeof(file_keymap_size_total));
            ifs.read((char *) &file_ctrl_block_byte_size, sizeof(file_ctrl_block_byte_size));
            ifs.read((char *) &file_keys_per_ctrl_block, sizeof(file_keys_per_ctrl_block));
            ifs.read((char *) &file_bitshift, sizeof(file_bitshift));
            ifs.read((char *) &file_values_size, sizeof(file_values_size));
            if (!ifs) InvalidIndex(name, "the file ends inside its header");

            size_t const expected_keymap_size = size_t{1} << m_main_bits;
            if (file_keymap_size != expected_keymap_size) {
                InvalidIndex(name, "key map of " + std::to_string(file_keymap_size) + " keys, this protal expects " +
                                   std::to_string(expected_keymap_size) + " (k = " + std::to_string(m_exact_k) + ")");
            }
            if (file_keys_per_ctrl_block == 0 || file_bitshift >= 64 || (size_t{1} << file_bitshift) != file_keys_per_ctrl_block) {
                InvalidIndex(name, "inconsistent control block layout (" + std::to_string(file_keys_per_ctrl_block) +
                                   " keys per block, shift " + std::to_string(file_bitshift) + ")");
            }
            if (file_ctrl_block_byte_size != sizeof(uint64_t)) {
                InvalidIndex(name, "control blocks of " + std::to_string(file_ctrl_block_byte_size) + " bytes, expected 8");
            }
            size_t const cell_size = file_ctrl_block_byte_size / sizeof(KeyMap_t);
            size_t const expected_total = file_keymap_size + ((file_keymap_size / file_keys_per_ctrl_block) + 1) * cell_size;
            if (file_keymap_size_total != expected_total) {
                InvalidIndex(name, "key map of " + std::to_string(file_keymap_size_total) + " cells, the layout implies " +
                                   std::to_string(expected_total));
            }
            // More values than the 64-bit address space holds is a corrupt header, not an allocation to try.
            if (file_values_size > (uint64_t{1} << 56) / sizeof(ValueEntry)) {
                InvalidIndex(name, std::to_string(file_values_size) + " values (corrupt file?)");
            }

            keymap_size = file_keymap_size;
            keymap_size_total = file_keymap_size_total;
            ctrl_block_byte_size = file_ctrl_block_byte_size;
            m_keys_per_ctrl_block = file_keys_per_ctrl_block;
            ctrl_block_frequency_bitshift = file_bitshift;
            values_size = file_values_size;
            // Fields derived from the ones above, recomputed so they match the loaded layout.
            ctrl_block_key_mask = (uint64_t{1} << ctrl_block_frequency_bitshift) - 1;
            ctrl_block_cell_size = cell_size;
            ctrl_block_byte_size_shift = log2(ctrl_block_cell_size);
            max_key_ubiquity = (1 << (sizeof(KeyMap_t)*8)) / m_keys_per_ctrl_block;
            return uint64_t{keymap_size_total} * sizeof(KeyMap_t) + uint64_t{values_size} * sizeof(ValueEntry);
        }


        void CountUpKey(uint64_t key) {
            assert(key < keymap_size);
            auto index = KeymapIndex(key);
            m_keymap[index] += m_keymap[index] < UINT16_MAX;
        }

        void CountUpKeyOMP(uint64_t key) {
            assert(key < keymap_size);
            auto index = KeymapIndex(key);

#pragma omp atomic
            m_keymap[index] += m_keymap[index] < UINT16_MAX;
        }

        uint64_t GetKey(uint64_t key) {
            return m_keymap[KeymapIndex(key)];
        }

        uint8_t Key(uint64_t key) {
            return m_keymap[KeymapIndex(key)];
        }

        // The 8-byte layout's values of a key (--build and the database tools; a query run packs its
        // values, GetPacked).
        void Get(uint64_t key, ValueEntry* &start, ValueEntry* &end, uint32_t* &flexblock_begin, uint32_t* &flexblock_end) {
            size_t main_key = MainKey(key);
            size_t flex_key = FlexKey(key);

            if (m_map == nullptr) {
                std::cerr << "The index is packed (Seedmap::Pack); this operation needs its 8-byte layout" << std::endl;
                exit(8);
            }
            if (main_key > m_keymask) {
                std::cout << "Key is larger than keymask" << key << " > " << m_keymask << std::endl;
            }
            uint64_t block_start_idx = ControlBlockIndex(main_key);
            uint64_t block_end_idx = block_start_idx + m_keys_per_ctrl_block + ctrl_block_cell_size;

            uint64_t block_value_start_idx = *((uint64_t *) (m_keymap + block_start_idx));
            uint64_t block_value_end_idx = *((uint64_t *) (m_keymap + block_end_idx));
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            size_t key_index = block_start_idx + ctrl_block_cell_size + BlockKey(main_key);

            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx +
                                                                                          m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            flexblock_begin = nullptr;
            flexblock_end = nullptr;
            if (key_value_block_size >= m_flex_threshold) {
                flexblock_begin = (uint32_t*) (m_map + key_value_start);
                flexblock_end = flexblock_begin + key_value_block_size - FlexBlockSize(key_value_block_size);
                key_value_start += FlexBlockSize(key_value_block_size);
            }

            start = m_map + key_value_start;
            end =  m_map + key_value_end;
        }

        void Get(uint64_t key, ValueEntry* &start, ValueEntry* &end) {
            size_t main_key = MainKey(key);
            size_t flex_key = FlexKey(key);

            if (main_key > m_keymask) {
                std::cout << "Key is larger than keymask" << key << " > " << m_keymask << std::endl;
            }
            uint64_t block_start_idx = ControlBlockIndex(main_key);
            uint64_t block_end_idx = block_start_idx + m_keys_per_ctrl_block + ctrl_block_cell_size;

            uint64_t block_value_start_idx = *((uint64_t *) (m_keymap + block_start_idx));
            uint64_t block_value_end_idx = *((uint64_t *) (m_keymap + block_end_idx));
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            size_t key_index = block_start_idx + ctrl_block_cell_size + BlockKey(main_key);

            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx +
                                                                                          m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                start = nullptr;
                end = nullptr;
                return;
            }

            if (key_value_block_size >= m_flex_threshold) {
                key_value_start += FlexBlockSize(key_value_block_size);
            }
//            if (key_value_block_size >= m_flex_threshold) {
//                std::cout << key_value_block_size << std::endl;
//                PrintBlock(main_key);
//
//                uint32_t* flexmer_ptr = (uint32_t*) (m_map + key_value_start);
//                // Update.. also give out positions of extended keys.
//                key_value_start += FlexBlockSize(key_value_block_size);
//
//                auto values_size = key_value_end - key_value_start;
//
//                std::cout << "Flexmer: " << KmerUtils::ToString(flex_key, m_flex_k_bits) << " " << SeedmapUtils::BitString<64>(flex_key) << std::endl;
//                for (auto i = 0; i < values_size; i++) {
//                    std::cout << " -> " << KmerUtils::ToString(flexmer_ptr[i], m_flex_k_bits) << " " << Similarity(flex_key, flexmer_ptr[i]) << " " << SeedmapUtils::BitString<64>(flexmer_ptr[i]) << std::endl;
//                }
//
//                Utils::Input();
//            }


            if (key_value_start > key_value_end) {
                std::cout << "Block end offset: " << m_keys_per_ctrl_block + ctrl_block_cell_size << std::endl;
                std::cout << "block_value_start_idx: " << block_value_start_idx << std::endl;
                std::cout << "block_value_end_idx: " << block_value_end_idx << std::endl;
                std::cout << "block_start_idx: " << block_start_idx << std::endl;
                std::cout << "key_index: " << key_index << std::endl;
                std::cout << "BlockKey(main_key): " << BlockKey(main_key) << std::endl;
                std::cout << main_key << " " << KmerUtils::ToString(main_key, m_main_bits) << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "key_value_start: " << key_value_start << std::endl;
                std::cout << "key_value_end: " << key_value_end << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
//                PrintBlock(main_key);
                std::cout << "Seedmap: 240" << std::endl;
                Print();
                exit(71);
            }

            start = m_map + key_value_start;
            end =  m_map + key_value_end;
        }

        bool Put(uint64_t &key, uint64_t &taxid, uint64_t &geneid, uint64_t &genepos) {
            if (m_map == nullptr) {
                std::cerr << "You need to initialize the values with BuildValuePointers()" << std::endl;
                exit(8);
            }
            if (key > keymap_max) {
                std::cerr << "Key too large (" << key << " > " << keymap_max << ")" << std::endl;
                exit(8);
            }


            uint64_t block_start_idx = ControlBlockIndex(key);
            uint64_t block_end_idx = block_start_idx + m_keys_per_ctrl_block + ctrl_block_cell_size;

            uint64_t block_value_start_idx = *((uint64_t*) (m_keymap + block_start_idx ) );
            uint64_t block_value_end_idx = *((uint64_t*) (m_keymap + block_end_idx) );
            uint64_t block_value_size = block_value_end_idx - block_value_start_idx;

            if (!block_value_size) {
                return false;
            }

            size_t key_index = block_start_idx + ctrl_block_cell_size + BlockKey(key);
            auto key_value_start = block_value_start_idx + m_keymap[key_index];
            auto key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx + m_keymap[key_index + 1];
            auto key_value_block_size = key_value_end - key_value_start;

            if (!key_value_block_size) {
                return false;
            }

            auto index = key_value_start;

            if (key_value_start > key_value_end) {
                std::cout << "block_value_start_idx: " << block_value_start_idx << std::endl;
                std::cout << "block_value_end_idx: " << block_value_end_idx << std::endl;
                std::cout << "block_start_idx: " << block_start_idx << std::endl;
                std::cout << "key_index: " << key_index << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;
                std::cout << "key_value_start: " << key_value_start << std::endl;
                std::cout << "key_value_end: " << key_value_end << std::endl;
                std::cout << "key_value_block_size: " << key_value_block_size << std::endl;


//                PrintBlock(key);
                std::cout << "Key value start (" << key_value_start << ") can not be larger than key value end (" << key_value_end << ")" << std::endl;
                exit(72);
            }

            for (; !m_map[index].Empty() && index < key_value_end; index++);

            if (index == key_value_end) {
                std::cout << "bucket too large, not present " << key_value_start << " size: " << (key_value_end - key_value_start) << std::endl;
                return false;
            }
            m_map[index].Put(taxid, geneid, genepos);

            ValueEntry *begin, *end;
            Get(key, begin, end);

            if constexpr(false) {
//            std::cout << "put index: " << index << " pointer: " << &m_map[index] << " begin pointer " << begin  << " " << taxid << " " << genepos << " " << geneid << std::endl;
                bool found = false;
                for (auto it = begin; it < end; it++) {
//                std::cout << it->ToString() << std::endl;
                    auto [t, g, p] = it->Get();
                    if (taxid == t && geneid == g && genepos == p) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cout << "____searched: " << taxid << ", " << geneid << ", " << genepos << std::endl;
                    for (auto it = begin; it < end; it++) {
                        std::cout << it->ToString() << std::endl;
                    }
                    exit(9);
                } else {
                    m_found_counter++;
                    if (m_found_counter % 10'000'000 == 0)
                        std::cout << "found counter: " << m_found_counter << std::endl;
                }
            }

//            std::string stop;
//            std::cin >> stop;
            return true;
        }

        uint64_t MainKey(uint64_t key) const {
            return (key & m_main_key_mask) >> m_flex_k ; // shift by half of flexbits which is flex_k
        }

        // Asks the CPU to fetch what Get(key, ...) reads of the key map: the key's control block and the next
        // block's value offset, 32 bytes from the block's start (one or two cache lines). The key map is gigabytes,
        // so a lookup misses the caches; fetched ahead, a read's lookups wait for memory together, not in turn.
        void PrefetchKey(uint64_t key) {
            char const* const block = reinterpret_cast<char const*>(m_keymap + ControlBlockIndex(MainKey(key)));
            __builtin_prefetch(block);
            __builtin_prefetch(block + 31);
        }

        uint64_t FlexKey(uint64_t key) const {
            return ((key & m_flex_key_mask_left) >> (m_main_bits)) | (key & m_flex_key_mask_right) ; // shift by half of flexbits which is flex_k
        }

        uint64_t FlexBlockSize(uint64_t block_size) const {
            return (block_size + 2) / 3;
        }

//        static uint32_t Similarity(uint32_t a, uint32_t b) {
//            return std::popcount(((~(a ^ b) >> 1) & ~(a ^ b)) & 0b00000000010101010101010101010101);
//        }
        static uint32_t Similarity(uint32_t a, uint32_t b) {
            return __builtin_popcount(((~(a ^ b) >> 1) & ~(a ^ b)) & 0b01010101010101010101010101010101);
//            return std::popcount(((~(a ^ b) >> 1) & ~(a ^ b)) & 0b01010101010101010101010101010101);
        }

        bool PutOMP(uint64_t &key, uint64_t &taxid, uint64_t &geneid, uint64_t &genepos) {
            ValueEntry entry;
            entry.Put(taxid, geneid, genepos + m_flex_k_half);  // the core sits flex_k/2 bases into the k-mer
            bool placed = false;
#pragma omp critical(put)
            placed = PutOwned(key, entry.value);
            return placed;
        }

        // Puts a value (the bits of a ValueEntry::Put) into the first empty slot of its key, and in a
        // flex block the key's flex part into the flex cell of that slot, in the 8-byte layout or the packed
        // one (BuildValuePointers with a layout, as --build fills them). No lock: one thread at a time
        // per control block, as PutOMP's lock or the index build's key ranges (Build.h) ensure; a
        // key's values arrive in reference order and fill its slots in that order. false for a key
        // without values (not counted, or dropped by BuildValuePointers).
        bool PutOwned(uint64_t key, uint64_t value) {
            uint64_t const main_key = MainKey(key);
            uint64_t const flex_key = FlexKey(key);
            if (m_map == nullptr && m_packed == nullptr) {
                std::cerr << "You need to initialize the values with BuildValuePointers()" << std::endl;
                exit(8);
            }
            if (main_key > keymap_max) {
                std::cerr << "Key too large (" << main_key << " > " << keymap_max << ")" << std::endl;
                exit(8);
            }

            uint64_t const block_start_idx = ControlBlockIndex(main_key);
            uint64_t const block_end_idx = block_start_idx + m_keys_per_ctrl_block + ctrl_block_cell_size;
            uint64_t const block_value_start_idx = *((uint64_t*) (m_keymap + block_start_idx));
            uint64_t const block_value_end_idx = *((uint64_t*) (m_keymap + block_end_idx));
            if (block_value_end_idx == block_value_start_idx) return false;

            size_t const key_index = block_start_idx + ctrl_block_cell_size + BlockKey(main_key);
            uint64_t const key_value_start = block_value_start_idx + m_keymap[key_index];
            uint64_t const key_value_end = (key_index + 1) == block_end_idx ? block_value_end_idx : block_value_start_idx + m_keymap[key_index + 1];
            uint64_t const key_value_block_size = key_value_end - key_value_start;
            if (!key_value_block_size) return false;

            bool const is_flex = key_value_block_size >= m_flex_threshold;
            uint64_t const flex_block_size = FlexBlockSize(key_value_block_size);
            uint64_t value_index = key_value_start + (is_flex * flex_block_size);
            auto failed = [&](char const* what, int code) {
                std::cerr << "Cannot place a value of key " << main_key << " (" << KmerUtils::ToString(main_key, m_main_bits)
                          << ", flex part " << KmerUtils::ToString(flex_key, m_flex_k_bits) << "): " << what << "; its values "
                          << key_value_start << "-" << key_value_end << (is_flex ? " (flex block of " : " (")
                          << (is_flex ? flex_block_size : 0) << "), " << values_size << " values in all" << std::endl;
                exit(code);
            };
            if (value_index >= values_size) failed("they lie beyond the values", 73);

            if (IsPacked()) {
                // The key's region of the packed values (PackedLayout): its flex cells, then its entries. They fill
                // in order, so the first empty entry follows the filled ones (found by bisection; a placed entry is
                // never 0, its unique flag is set). Bytes shared with the neighbouring keys, which another thread
                // may be filling, are read and written atomically (GetBits, PutBits).
                uint64_t const region = key_value_start * m_slot_bits, region_end = key_value_end * m_slot_bits;
                uint64_t const safe_begin = (region >> 3) + ((region & 7) ? 1 : 0), safe_end = region_end >> 3;
                uint64_t const entries = key_value_block_size - (is_flex ? flex_block_size : 0);
                uint64_t const entry_bit = region + (is_flex ? 32 * entries : 0);
                uint64_t lo = 0, hi = entries;
                while (lo < hi) {
                    uint64_t const mid = lo + (hi - lo) / 2;
                    if (GetBits(m_packed, entry_bit + mid * m_entry_bits, m_entry_bits, safe_begin, safe_end) != 0) lo = mid + 1;
                    else hi = mid;
                }
                if (lo == entries) failed("all its slots are taken", 75);
                uint64_t packed = 0;
                if (!PackValue(value, packed)) failed("its taxid, gene id or position is wider than the reference's (PackedLayout)", 8);
                PutBits(m_packed, entry_bit + lo * m_entry_bits, packed, m_entry_bits, safe_begin, safe_end);
                if (is_flex) {
                    uint64_t const flex_bit = region + 32 * lo;
                    if (GetBits(m_packed, flex_bit, 32, safe_begin, safe_end) != 0) failed("its flex cell is taken", 76);
                    PutBits(m_packed, flex_bit, static_cast<uint32_t>(flex_key), 32, safe_begin, safe_end);
                }
                return true;
            }

            // The next empty slot
            while (value_index < key_value_end && !m_map[value_index].Empty()) value_index++;
            if (value_index == key_value_end) failed("all its slots are taken", 75);
            m_map[value_index].value = value;

            if (is_flex) {
                uint64_t const flexpos = value_index - key_value_start - flex_block_size;
                uint32_t* const flex_ptr = (uint32_t*) &(m_map[key_value_start]);
                if (flex_ptr[flexpos] != 0) failed("its flex cell is taken", 76);
                flex_ptr[flexpos] = static_cast<uint32_t>(flex_key);
            }
            return true;
        }


        struct subkey {
            uint32_t key = 0;
            uint32_t count = 0;

            void Reset() {
                key = 0; count = 0;
            }
            void Print() {
                std::cout << "key: " << key << "  count: " << count << std::endl;
            }
            std::string ToString() const {
                return std::to_string(key) + '\t' + std::to_string(count);
            }
        };

        void BuildValuePointersFirstIteration(uint64_t& total_count, subkey* subkey, int* kmer_frequencies,
                                              uint64_t block_start_index) {
            // Set ubiquity of each key in subkey
            for (auto key_pos = 0; key_pos < m_keys_per_ctrl_block; key_pos++) {
                subkey[key_pos].key = key_pos;
                subkey[key_pos].count = m_keymap[block_start_index + key_pos];
                kmer_frequencies[subkey[key_pos].count]++;
                total_count += subkey[key_pos].count;
            }

            // Sort so keys are sorted by ubiquity.
            std::sort(subkey, subkey + m_keys_per_ctrl_block, [](const struct subkey &a, const struct subkey &b) {
                return a.count < b.count;
            });
        }

        void BuildValuePointersSecondIteration(subkey* subkey, size_t& count_failed_demand, size_t& count_total_stored,
                                               size_t& count_total_demand, uint64_t& max_keyblock_size) {
            uint64_t current_sum = 0;
            uint64_t total_demand = 0;

            for (auto key_pos = 0; key_pos < m_keys_per_ctrl_block; key_pos++) {
                if (subkey[key_pos].count == 0) continue;

                // Given the Flex-K approach block demand quantifies the number of value cells
                // a given key needs.
                size_t key_value_space_demand = subkey[key_pos].count >= m_flex_threshold ?
                                                subkey[key_pos].count + DivisionByTwoCeiling(subkey[key_pos].count) : subkey[key_pos].count;

                total_demand += key_value_space_demand;
//                    subkey[key_pos].Print();
                // This if decides whether keys are taken or not !!!
                // Important to get the key_value_space demand correct otherwise the offsets will be wrong
                // and it will come to SEGMENTATION FAULTS
                if (subkey[key_pos].count < max_key_multiplicity && key_value_space_demand < max_keyblock_size && (current_sum + key_value_space_demand) < max_block_size) {
                    current_sum += key_value_space_demand;
                } else {
                    subkey[key_pos].count = 0;
                }
            }
            // Only statistics.
            count_total_demand += total_demand;
            count_total_stored += current_sum;
            count_failed_demand += (current_sum < total_demand);

            std::sort(subkey, subkey + m_keys_per_ctrl_block, [](const struct subkey &a, const struct subkey &b) {
                return a.key < b.key;
            });
        }

        void BuildValuePointersBlockThirdIteration(subkey* subkey, uint64_t block_start_index, uint64_t& global_position) {
            uint64_t local_position = 0;
            for (auto key_pos = 0; key_pos < m_keys_per_ctrl_block; key_pos++) {
                m_keymap[block_start_index + key_pos] = local_position;
                size_t key_value_space_demand = subkey[key_pos].count >= m_flex_threshold ?
                                                subkey[key_pos].count + DivisionByTwoCeiling(subkey[key_pos].count) : subkey[key_pos].count;

                local_position += key_value_space_demand;
                global_position += key_value_space_demand;
            }
        }

        void BuildValuePointersBlock(uint64_t block, uint64_t& global_position,
                                     size_t& count_failed_demand, size_t& count_total_stored, size_t& count_total_demand,
                                     int* kmer_frequencies, uint64_t& max_keyblock_size, subkey* subkey) {

            uint64_t* control_ptr = nullptr;
            auto block_index = block * (m_keys_per_ctrl_block + ctrl_block_cell_size);
            auto block_start_index = block_index + ctrl_block_cell_size;

            uint64_t total_count = 0;

            control_ptr = (uint64_t*) (m_keymap + block_index);
            *control_ptr = global_position;


            //####################################################################################################
            // first iteration
            // - After this keys in block are sorted by their frequency
            // second iteration
            // - After this, keys in block are resorted by their key.
            // Third iteration
            // - Here, the global offset values and local offset values are set.
            BuildValuePointersFirstIteration(total_count, subkey, kmer_frequencies, block_start_index);
            if (!total_count) return;
            BuildValuePointersSecondIteration(subkey, count_failed_demand, count_total_stored, count_total_demand,
                                              max_keyblock_size);
            BuildValuePointersBlockThirdIteration(subkey, block_start_index, global_position);
        }



        // A row of unique_kmers.tsv for each gene id of each reference taxon: row first_row[taxid] +
        // gene id - 1, so taxon t has first_row[t + 1] - first_row[t] gene ids (see Build.h).
        struct GeneRows {
            std::vector<uint64_t> first_row;
            uint64_t Size() const { return first_row.empty() ? 0 : first_row.back(); }
        };

        struct UniqueKmerTotals {
            size_t short_unique = 0, long_unique = 0, long_unique_two = 0, non_unique = 0;
            size_t comparisons = 0;  // flex parts compared for the distance-two flag
        };

        // Whether each of the n flex parts has another within distance 1 (all of the flex_k positions
        // equal but at most one). Two such parts agree once that position is masked out, so sorting the
        // parts with each position masked in turn finds them all: n log n per position instead of n^2.
        static void FlexNeighbours(uint32_t const* flex, size_t n, size_t flex_k, std::vector<uint8_t>& close,
                                   std::vector<std::pair<uint32_t, uint32_t>>& keyed) {
            close.assign(n, 0);
            keyed.resize(n);
            for (size_t p = 0; p < flex_k; p++) {
                uint32_t const mask = ~(3u << (2 * p));
                for (size_t e = 0; e < n; e++) keyed[e] = { flex[e] & mask, static_cast<uint32_t>(e) };
                std::sort(keyed.begin(), keyed.end());
                for (size_t a = 0, b; a < n; a = b) {
                    for (b = a + 1; b < n && keyed[b].first == keyed[a].first; b++) {}
                    if (b - a > 1) for (size_t e = a; e < b; e++) close[keyed[e].second] = 1;
                }
            }
        }

        // unique_kmers.tsv: for each gene with values in the index, its unique k-mers (a core with one
        // value: short; with several, unique by the flex part: long; long and at flex distance two or
        // more from every other value of the core: long, two) and all its k-mers, as counts and
        // fractions. Also sets the distance-two flag of those long unique values. Threads take ranges
        // of control blocks and count into one row per gene (rows); the rows are written in the order
        // of their genes' first values in the index, through a sparse_map filled in that order, which
        // is the order the serial scan with four string-keyed maps wrote them in. So the table does not
        // depend on the threads, and is the one earlier builds wrote.
        UniqueKmerTotals CountUniqueKmers(std::ostream& os, GeneRows const& rows, int threads) {
            struct Counts { uint32_t short_unique = 0, long_unique = 0, long_unique_two = 0, total = 0; };
            uint64_t const n_rows = rows.Size();
            std::vector<Counts> counts(n_rows);
            std::vector<uint64_t> first_value(n_rows, UINT64_MAX);
            UniqueKmerTotals totals;
            size_t outside = 0;
            int64_t const n_blocks = static_cast<int64_t>((keymap_max + m_keys_per_ctrl_block - 1) / m_keys_per_ctrl_block);
            // From this many values of a core on, FlexNeighbours is cheaper than comparing pairs.
            constexpr size_t kSortFrom = 256;

#pragma omp parallel num_threads(std::max(threads, 1))
            {
                UniqueKmerTotals local;
                size_t local_outside = 0;
                std::vector<uint8_t> close;
                std::vector<std::pair<uint32_t, uint32_t>> keyed;
                std::vector<uint32_t> flex_cells;  // a packed key's flex cells
                PackedBlock packed;

#pragma omp for schedule(dynamic, 1 << 14)
                for (int64_t block = 0; block < n_blocks; block++) {
                    uint64_t const ctrl = ControlBlockIndex(static_cast<uint64_t>(block) * m_keys_per_ctrl_block);
                    uint64_t const begin = *((uint64_t*)(m_keymap + ctrl));
                    uint64_t const end = *((uint64_t*)(m_keymap + ctrl + m_keys_per_ctrl_block + ctrl_block_cell_size));
                    if (begin == end) continue;

                    for (size_t j = ctrl_block_cell_size; j < ctrl_block_cell_size + m_keys_per_ctrl_block; j++) {
                        uint64_t const key_start = begin + m_keymap[ctrl + j];
                        uint64_t const key_end = j + 1 == ctrl_block_cell_size + m_keys_per_ctrl_block ? end : begin + m_keymap[ctrl + j + 1];
                        uint64_t const size = key_end - key_start;
                        if (size == 0) continue;
                        bool const has_flex = size >= m_flex_threshold;
                        uint64_t const flex_slots = has_flex ? FlexBlockSize(size) : 0;
                        uint64_t const n = size - flex_slots;
                        // The key's values in either layout (--build holds them packed).
                        ValueEntry* const values = IsPacked() ? nullptr : m_map + key_start + flex_slots;
                        if (IsPacked()) BlockOfSlots(key_start, size, packed);
                        auto entry_of = [&](size_t e) {
                            ValueEntry entry;
                            entry.value = values ? values[e].value : EntryValue(packed, static_cast<uint32_t>(e));
                            return entry;
                        };
                        auto set_two = [&](size_t e) {
                            if (values) values[e].SetFlagUniqueDistanceMinTwo();
                            else SetUniqueDistanceTwoFlag(packed, static_cast<uint32_t>(e));
                        };

                        if (has_flex) {
                            // Flex part e belongs to value e. Only unique values get the flag, and only
                            // whether another value is within distance 1 matters.
                            uint32_t const* flex = values ? (uint32_t const*)(m_map + key_start) : nullptr;
                            if (!values) {
                                flex_cells.resize(n);
                                for (size_t e = 0; e < n; e++) flex_cells[e] = FlexCell(packed, static_cast<uint32_t>(e));
                                flex = flex_cells.data();
                            }
                            if (n >= kSortFrom) {
                                FlexNeighbours(flex, n, m_flex_k, close, keyed);
                                local.comparisons += n * m_flex_k;
                                for (size_t e = 0; e < n; e++) {
                                    if (entry_of(e).IsFlagUnique() && !close[e]) set_two(e);
                                }
                            } else {
                                for (size_t e = 0; e < n; e++) {
                                    if (!entry_of(e).IsFlagUnique()) continue;
                                    bool near = false;
                                    for (size_t o = 0; o < n && !near; o++) {
                                        if (o == e) continue;
                                        local.comparisons++;
                                        near = Similarity(flex[e], flex[o]) + 1 >= m_flex_k;
                                    }
                                    if (!near) set_two(e);
                                }
                            }
                        }

                        for (size_t e = 0; e < n; e++) {
                            ValueEntry const entry = entry_of(e);
                            auto [taxid, geneid, pos] = entry.Get();
                            if (taxid + 1 >= rows.first_row.size() || geneid == 0 ||
                                geneid > rows.first_row[taxid + 1] - rows.first_row[taxid]) {
                                local_outside++;
                                continue;
                            }
                            uint64_t const row = rows.first_row[taxid] + geneid - 1;
                            bool const unique = entry.IsFlagUnique();
                            bool const two = unique && has_flex && entry.IsFlagUniqueDistanceMinTwo();
                            Counts& c = counts[row];
                            std::atomic_ref<uint32_t>(c.total).fetch_add(1, std::memory_order_relaxed);
                            if (unique) std::atomic_ref<uint32_t>(has_flex ? c.long_unique : c.short_unique).fetch_add(1, std::memory_order_relaxed);
                            if (two) std::atomic_ref<uint32_t>(c.long_unique_two).fetch_add(1, std::memory_order_relaxed);
                            uint64_t const position = key_start + flex_slots + e;
                            std::atomic_ref<uint64_t> first(first_value[row]);
                            uint64_t seen = first.load(std::memory_order_relaxed);
                            while (position < seen && !first.compare_exchange_weak(seen, position, std::memory_order_relaxed)) {}

                            local.short_unique += unique && !has_flex;
                            local.long_unique += unique && has_flex;
                            local.long_unique_two += two;
                            local.non_unique += !unique;
                        }
                    }
                }

#pragma omp critical(unique_kmer_totals)
                {
                    totals.short_unique += local.short_unique;
                    totals.long_unique += local.long_unique;
                    totals.long_unique_two += local.long_unique_two;
                    totals.non_unique += local.non_unique;
                    totals.comparisons += local.comparisons;
                    outside += local_outside;
                }
            }
            if (outside) {
                std::cerr << outside << " index values name a gene that is not in reference.map; rebuild the index" << std::endl;
                exit(8);
            }

            std::vector<uint64_t> order;
            for (uint64_t row = 0; row < n_rows; row++) {
                if (counts[row].total) order.push_back(row);
            }
            std::sort(order.begin(), order.end(), [&](uint64_t a, uint64_t b) { return first_value[a] < first_value[b]; });
            tsl::sparse_map<std::string, uint32_t> table;  // filled in that order; its iteration order is the file's
            for (auto row : order) {
                auto const taxid = static_cast<uint64_t>(
                        std::upper_bound(rows.first_row.begin(), rows.first_row.end(), row) - rows.first_row.begin()) - 1;
                table.insert({ std::to_string(taxid) + '_' + std::to_string(row - rows.first_row[taxid] + 1), 0 });
            }
            for (auto const& [key, _] : table) {
                auto [taxid, geneid] = KmerUtils::ExtractHeaderInformation(key);
                auto const& c = counts[rows.first_row[taxid] + geneid - 1];
                os << taxid << '\t'; //1
                os << geneid << '\t'; //2
                os << c.short_unique << '\t'; //3
                os << static_cast<double>(c.short_unique)/c.total << '\t'; //4
                os << c.long_unique << '\t'; //5
                os << static_cast<double>(c.long_unique)/c.total << '\t'; //6
                os << c.long_unique_two << '\t';
                os << static_cast<double>(c.long_unique_two)/c.total << '\t';
                os << c.total << '\n';
            }
            std::cout << "ShortUniques:    " << totals.short_unique << std::endl;
            std::cout << "LongUniques:    " << totals.long_unique << std::endl;
            std::cout << "LongUniquesTwo:    " << totals.long_unique_two << std::endl;
            std::cout << "Nonuniques: " << totals.non_unique << std::endl;
            return totals;
        }

        // The value positions of every control block and key (from the counts of pass 1), in
        // `threads` threads: each part of the blocks is laid out from 0, then the parts' sizes are
        // summed in order and each part's start is added to its blocks' positions. The same key
        // map as on one thread. With a layout the values are allocated packed (PackedLayout), as a query
        // run holds them, and filled so (PutOwned): --build never holds the 8-byte layout (36 GB at r226
        // instead of 29 with the key map), which only the files keep.
        void BuildValuePointers(int threads = 1, PackedLayout const* pack = nullptr) {
            uint64_t const num_ctrl_blocks = keymap_size >> ctrl_block_frequency_bitshift;
            std::cout << "number ctrl blocks: " << num_ctrl_blocks << std::endl;
            uint64_t const max_keyblock_size = max_key_ubiquity*2;
            uint64_t const block_cells = m_keys_per_ctrl_block + ctrl_block_cell_size;

            // Understand kmer frequencies better to improve sensitivity of alignment.
            constexpr int kmer_freq_size = 1 << (sizeof(KeyMap_t) * 8);
            std::vector<int> kmer_frequencies(kmer_freq_size, 0);
            size_t count_failed_demand = 0;
            size_t count_total_stored = 0;
            size_t count_total_demand = 0;

            threads = std::max(threads, 1);
            size_t const parts = threads == 1 ? 1 : static_cast<size_t>(threads) * 16;
            std::vector<uint64_t> part_start(parts + 1, 0);
            auto part_begin = [&](size_t part) { return num_ctrl_blocks * part / parts; };

#pragma omp parallel num_threads(threads)
            {
                std::vector<int> frequencies(kmer_freq_size, 0);
                size_t failed_demand = 0, total_stored = 0, total_demand = 0;
                uint64_t keyblock_size = max_keyblock_size;
                subkey keys[8];
                std::vector<subkey> key_storage;
                subkey* sk = keys;
                if (m_keys_per_ctrl_block > 8) {
                    key_storage.resize(m_keys_per_ctrl_block);
                    sk = key_storage.data();
                }

                // Iterate through blocks. Each block manages <m_keys_per_ctrl_block> 15-mers
#pragma omp for schedule(dynamic, 1)
                for (size_t part = 0; part < parts; part++) {
                    uint64_t position = 0;
                    for (uint64_t block = part_begin(part); block < part_begin(part + 1); block++) {
                        BuildValuePointersBlock(block, position, failed_demand, total_stored, total_demand,
                                                frequencies.data(), keyblock_size, sk);
                    }
                    part_start[part + 1] = position;
                }

#pragma omp critical(value_pointers)
                {
                    for (int i = 0; i < kmer_freq_size; i++) kmer_frequencies[i] += frequencies[i];
                    count_failed_demand += failed_demand;
                    count_total_stored += total_stored;
                    count_total_demand += total_demand;
                }
            }
            for (size_t part = 0; part < parts; part++) part_start[part + 1] += part_start[part];
            if (parts > 1) {
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
                for (size_t part = 1; part < parts; part++) {
                    for (uint64_t block = part_begin(part); block < part_begin(part + 1); block++) {
                        *((uint64_t*) (m_keymap + block * block_cells)) += part_start[part];
                    }
                }
            }
            uint64_t const global_position = part_start[parts];

            // set last control pointer
            uint64_t* control_ptr = (uint64_t*) (m_keymap + keymap_size_total - ctrl_block_cell_size);
            *control_ptr = global_position;
            values_size = global_position;
            std::cout << "ctrl_block_cell_size: " << ctrl_block_cell_size << std::endl;

            if (pack) {
                SetLayout(*pack);
                AllocatePacked();  // zeroed: an all-zero entry is empty
                uint64_t entries = 0;
#pragma omp parallel for num_threads(threads) schedule(dynamic, 1 << 14) reduction(+:entries)
                for (int64_t block = 0; block < static_cast<int64_t>(num_ctrl_blocks); block++) {
                    uint64_t const begin = BlockStart(block), end = BlockStart(block + 1);
                    if (begin == end) continue;
                    for (size_t j = 0; j < m_keys_per_ctrl_block; j++) {
                        uint64_t const start = begin + m_keymap[block * block_cells + ctrl_block_cell_size + j];
                        uint64_t const stop = j + 1 == m_keys_per_ctrl_block ? end : begin + m_keymap[block * block_cells + ctrl_block_cell_size + j + 1];
                        uint64_t const S = stop - start;
                        entries += S >= m_flex_threshold ? S - FlexBlockSize(S) : S;
                    }
                }
                m_packed_entries = entries;
            } else {
                AllocateValues(values_size, true);  // zeroed: an all-zero entry is empty
            }
            std::cout << "Index in memory: " << MemoryDescription() << std::endl;

            constexpr bool verbose = true;
            if constexpr(verbose) {
                std::cout << "Values size: " << values_size << std::endl;
                std::cout << "space requirements" << std::endl;
                double key_mem = (double) keymap_size_total / (1024 * 1024 * 1024);
                double val_mem = (double) (values_size * 8) / (1024 * 1024 * 1024);
                std::cout << "keys:   " << key_mem << " GB" << std::endl;
                std::cout << "values: " << val_mem << " GB" << std::endl;
                std::cout << "total:  " << key_mem + val_mem << " GB" << std::endl;


                for (auto i = 0; i < kmer_freq_size; i++) {
                    std::cout << i << '\t' << kmer_frequencies[i] << '\n';
                }
                std::cout << "Total stored: " << count_total_stored << std::endl;
                std::cout << "Total demand: " << count_total_demand << std::endl;
                std::cout << "Stored "
                          << (100 * static_cast<double>(count_total_stored) / static_cast<double>(count_total_demand))
                          << " of total values." << std::endl;
                std::cout << "Failed to meet demands? " << count_failed_demand << std::endl;
            }
        }


        void SortForKeys() {
            std::cout << "Sort for keys" << std::endl;
            size_t row_index = 0;
            std::string line = "";
            size_t non_empty_block_count = 0;
            for (auto key = 0; key < keymap_max; key++) {
                auto ctrl_block = ControlBlockIndex(key);
                auto value_block_idx_start = *((uint64_t*) (m_keymap + ctrl_block));
                auto value_block_idx_end = *((uint64_t*) (m_keymap + ctrl_block + m_keys_per_ctrl_block + ctrl_block_cell_size));
                auto value_block_size = value_block_idx_end - value_block_idx_start;

                std::sort(m_map+value_block_idx_start,m_map+value_block_idx_end);
            }
        }

// ctrl_block_keys_index
//  │                       Values array
//  │   │          │       │             │
//  │   ├──────────┤       ├─────────────┤
//  └──►│CTRL-BLOCK├──────►│             │<- ctrl_block_values_begin
//      ├──────────┤       │             │
//      │Key1      │       │             │
//      ├──────────┤       │             │
//      │Key2      │       │             │
//      ├──────────┤       │             │
//      │Key3      │       │             │
//      ├──────────┤       │             │    m_keys_per_ctrl_block
//      │Key4      │       │             │
//      ├──────────┤       │             │    #keys in this range
//      │Key5      │       │             │
//      ├──────────┤       │             │
//      │Key6      │       │             │
//      ├──────────┤       │             │
//      │Key7      │       │             │
//      ├──────────┤       │             │
//      │Key8      │       │             │
//      ├──────────┤       │             │
//      │CTRL-Block├───┐   │             │
//      ├──────────┤   │   │             │
//      │          │   │   ├─────────────┤
//      │          │   └──►│             │<- ctrl_block_values_end
//      ├──────────┤       │             │
//      │          │

        void PrintBlock(size_t key) {
            std::cout << "Print block key: " << KmerUtils::ToString(key, m_main_bits) << std::endl;
            auto ctrl_block_keys_index = ControlBlockIndex(key);
            auto ctrl_block_values_begin = *((uint64_t*)(m_keymap + ctrl_block_keys_index));
            auto ctrl_block_values_end = *((uint64_t*)(m_keymap + ctrl_block_keys_index + m_keys_per_ctrl_block + ctrl_block_cell_size));
            auto ctrl_block_values_size = ctrl_block_values_end - ctrl_block_values_begin;

            std::cout << "\n#######################################\n## KEY ARRAY ###########################~\n" << std::endl;
            std::cout << "-- Ctrl Block --- Keys: " << ctrl_block_keys_index << " -- Values: " << ctrl_block_values_begin << std::endl;
            std::cout << "Idx in values: " << ctrl_block_values_begin << " - " << ctrl_block_values_end << std::endl;
            std::cout << "----------------" << std::endl;
            for (int i = ctrl_block_cell_size; i < ctrl_block_cell_size + m_keys_per_ctrl_block; i++) {
                std::cout << (i-ctrl_block_cell_size) << ": " << (uint32_t) m_keymap[ctrl_block_keys_index + i] << ", ";
            }
            std::cout << std::endl;

            std::cout << "-- Next Ctrl Block --- Keys: " << ctrl_block_keys_index + m_keys_per_ctrl_block + ctrl_block_cell_size << " -- Values: " << ctrl_block_values_end << std::endl;
            std::cout << "   .... " << std::endl;
            std::cout << "\n## VALUE ARRAY ###########################~" << std::endl;

            Utils::Input();

            for (int j = ctrl_block_cell_size; j < ctrl_block_cell_size + m_keys_per_ctrl_block; j++) {
                size_t key_index = ctrl_block_keys_index + j;
                auto key_value_start = ctrl_block_values_begin + m_keymap[key_index];
                auto key_value_end = j == ctrl_block_cell_size + m_keys_per_ctrl_block - 1 ?
                        ctrl_block_values_end : ctrl_block_values_begin + m_keymap[key_index + 1];
                auto key_value_block_size = key_value_end - key_value_start;
                auto i = 0;
                std::cout << "-- Value index: " << key_value_start << " -------------" << (j - ctrl_block_cell_size) << " From, To: " << key_value_start << " - " << key_value_end << " (Size: " << FlexBlockSize(key_value_block_size) << ")      ";
                if (key_value_block_size >= m_flex_threshold) {
                    std::cout << " Flexi-K block (" << key_value_start << ")" << std::endl;
                    auto flex_block_size = FlexBlockSize(key_value_block_size);
                    for (; i < flex_block_size; i++) {
                        auto* cell = m_map + key_value_start + i;
                        std::cout << 2*i << ": " << KmerUtils::ToString(((uint32_t*)cell)[0], m_flex_k_bits) << "  " << (2*i + 1) << ": " << KmerUtils::ToString(((uint32_t*)cell)[1], m_flex_k_bits) << std::endl;
                    }
                    std::cout << "-------------";
                }
                std::cout << " Value block" << std::endl;
                for (; i < key_value_block_size; i++) {
                    std::cout << key_value_start + i << ": " << m_map[key_value_start + i].ToString() << std::endl;
                }
            }
            std::cout << "-------------";
            std::cout << "End printing block \n##########################" << std::endl;
        }

        void PrintKeyMap(size_t print_n, size_t wrap_around, bool skip_empty=true) {
            size_t row_index = 0;
            std::string line = "";
            size_t non_empty_block_count = 0;
            for (auto row = 0; row_index < keymap_size_total; row++, row_index = row * wrap_around) {
                bool line_empty = true;
                // check if there is a non empty cell
                for (auto col = 0; col < wrap_around && row_index + col < keymap_size_total; col++) {
                    if (m_keymap[row_index + col]) {
//                        std::cout << "line not empty(" << row + col << ", " << row << ", " << col << ") " << (uint32_t) keymap[row_index + col] << std::endl;
                        line_empty = false;
                        break;
                    }
                }
                // print if there is non-empty cell
                if (line_empty) continue;
                non_empty_block_count++;
                if (non_empty_block_count > print_n) return;

                std::cout << row_index << " ";
                for (auto col = 0; col < wrap_around; col++) {
                    std::cout << SeedmapUtils::BitString<8>(m_keymap[row_index + col]) << " ";
                }
                std::cout << '\n';
            }
        }
    };
}
