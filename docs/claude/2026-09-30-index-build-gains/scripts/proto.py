"""Prototype (scratch only), on top of instrument.py: CountUniqueKmers counts per gene in a vector
found through an integer-keyed robin_map (not four string-keyed sparse_maps updated per value); the
string key goes once per gene into a sparse_map in the same order as before, so the table is written
in the same row order; the unused post-build Check pass is not run. Usage: proto.py <source dir>"""
import sys

src = sys.argv[1]


def edit(path, pairs):
    path = f"{src}/{path}"
    text = open(path).read()
    for old, new in pairs:
        if text.count(old) != 1:
            sys.exit(f"{path}: {old[:70]!r} found {text.count(old)} times")
        text = text.replace(old, new)
    open(path, "w").write(text)


edit("src/Hash/Seedmap.h", [
    ('#include "sparse_map.h"\n', '#include "sparse_map.h"\n#include "robin_map.h"\n'),
    ("            tsl::sparse_map<std::string, uint32_t> short_unique_kmers;\n"
     "            tsl::sparse_map<std::string, uint32_t> long_unique_kmers;\n"
     "            tsl::sparse_map<std::string, uint32_t> long_unique_two_kmers;\n"
     "            tsl::sparse_map<std::string, uint32_t> all_kmers;\n",
     "            struct GeneCounts { uint32_t short_u = 0, long_u = 0, long_two = 0, total = 0; };\n"
     "            tsl::sparse_map<std::string, uint32_t> row_order;  // gene -> row, inserted as before\n"
     "            tsl::robin_map<uint64_t, uint32_t> row_of;\n"
     "            std::vector<GeneCounts> gene_counts;\n"),
    ("                        std::string key_str = std::to_string(taxid) + '_' + std::to_string(geneid);\n",
     "                        uint64_t const gene_key = (uint64_t(taxid) << 32) | geneid;\n"
     "                        auto found = row_of.find(gene_key);\n"
     "                        uint32_t row;\n"
     "                        if (found == row_of.end()) {\n"
     "                            row = gene_counts.size();\n"
     "                            row_of.emplace(gene_key, row);\n"
     "                            row_order.insert({std::to_string(taxid) + '_' + std::to_string(geneid), row});\n"
     "                            gene_counts.emplace_back();\n"
     "                        } else {\n"
     "                            row = found->second;\n"
     "                        }\n"
     "                        auto& counts = gene_counts[row];\n"),
    ("                        if (!short_unique_kmers.contains(key_str)) {\n"
     "                            short_unique_kmers.insert({key_str, 0});\n"
     "                            long_unique_kmers.insert({key_str, 0});\n"
     "                            long_unique_two_kmers.insert({key_str, 0});\n"
     "                            all_kmers.insert({key_str, 0});\n"
     "                        }\n"
     "                        long_unique_two_kmers[key_str] += entry.IsFlagUnique() && entry.IsFlagUniqueDistanceMinTwo() && has_flex_block;\n"
     "                        long_unique_kmers[key_str] += entry.IsFlagUnique() && has_flex_block;\n"
     "                        short_unique_kmers[key_str] += entry.IsFlagUnique() && !has_flex_block;\n"
     "                        all_kmers[key_str]++;\n",
     "                        counts.long_two += entry.IsFlagUnique() && entry.IsFlagUniqueDistanceMinTwo() && has_flex_block;\n"
     "                        counts.long_u += entry.IsFlagUnique() && has_flex_block;\n"
     "                        counts.short_u += entry.IsFlagUnique() && !has_flex_block;\n"
     "                        counts.total++;\n"),
    ("            for (auto [key, short_uniques] : short_unique_kmers) {\n"
     "                auto long_uniques_two = long_unique_two_kmers[key];\n"
     "                auto long_uniques = long_unique_kmers[key];\n"
     "                auto total = all_kmers[key];\n",
     "            for (auto const& [key, row] : row_order) {\n"
     "                auto const& c = gene_counts[row];\n"
     "                auto short_uniques = c.short_u;\n"
     "                auto long_uniques_two = c.long_two;\n"
     "                auto long_uniques = c.long_u;\n"
     "                auto total = c.total;\n"),
    ('<< " genes=" << short_unique_kmers.size()', '<< " genes=" << gene_counts.size()'),
    # all pairs: only unique entries read their distance, and only whether it is above 1
    ("                        for (uint32_t* flex_cell = flex_begin; flex_cell != flex_end; flex_cell++) {\n"
     "                            auto max = 0;\n"
     "                            for (uint32_t* flex_cell_sim = flex_begin; flex_cell_sim != flex_end; flex_cell_sim++) {\n"
     "                                auto sim = Seedmap::Similarity(*flex_cell, *flex_cell_sim);\n"
     "                                max = sim > max && flex_cell != flex_cell_sim ? sim : max;\n"
     "                            }\n",
     "                        auto* block_entries = m_map + key_value_start + FlexBlockSize(key_value_block_size);\n"
     "                        for (uint32_t* flex_cell = flex_begin; flex_cell != flex_end; flex_cell++) {\n"
     "                            uint32_t max = 0;\n"
     "                            if (block_entries[flex_cell - flex_begin].IsFlagUnique()) {\n"
     "                                for (uint32_t* other = flex_begin; other != flex_end && max + 1 < m_flex_k; other++) {\n"
     "                                    if (other == flex_cell) continue;\n"
     "                                    auto sim = Seedmap::Similarity(*flex_cell, *other);\n"
     "                                    max = sim > max ? sim : max;\n"
     "                                }\n"
     "                            }\n"),
])

edit("src/Hash/KmerLookup.h", [
    ("        inline void GetFlex(size_t &kmer, std::vector<ValueEntry*>& max_sim_entries, uint32_t& max_similarity) {\n",
     "        // The entries whose whole k-mer is kmer (none if more than m_max_ubiquity are): what GetFlex\n"
     "        // gives when the best similarity is the whole flex part, the only case the uniqueness check uses.\n"
     "        inline void GetExact(size_t &kmer, std::vector<ValueEntry*>& entries) {\n"
     "            m_sm.Get(kmer, m_entry_begin, m_entry_end, m_flex_begin, m_flex_end);\n"
     "            I_queries++; if (m_entry_begin != nullptr && m_flex_begin != nullptr) I_scanned += m_flex_end - m_flex_begin;\n"
     "            else if (m_entry_begin != nullptr && m_entry_end - m_entry_begin == 1) I_single++;\n"
     "            if (m_entry_begin == nullptr || m_entry_end == nullptr || m_flex_begin == nullptr) return;\n"
     "            uint32_t const flex_key = static_cast<uint32_t>(m_sm.FlexKey(kmer));\n"
     "            for (auto cell = m_flex_begin; cell < m_flex_end; cell++) {\n"
     "                if (*cell == flex_key) entries.emplace_back(m_entry_begin + (cell - m_flex_begin));\n"
     "            }\n"
     "            if (entries.size() > m_max_ubiquity) entries.clear();\n"
     "        }\n\n"
     "        inline void GetFlex(size_t &kmer, std::vector<ValueEntry*>& max_sim_entries, uint32_t& max_similarity) {\n"),
])

edit("src/Build.h", [
    ("                lookup.GetFlex(pair.first, max_sim_entries, max_sim);\n                if (max_sim_entries.empty()) {\n",
     "                lookup.GetExact(pair.first, max_sim_entries);\n"
     "                max_sim = max_sim_entries.empty() ? 0 : best_possible_sim;\n"
     "                if (max_sim_entries.empty()) {\n"),
])

edit("src/RunProtal.h", [
    ("            protal::build::Check<SimpleKmerHandler<ClosedSyncmer>, KmerPutterSM, DEBUG_NONE>(\n"
     "                    options, kmer_putter, iterator);\n", ""),
])
print("prototype", src)
