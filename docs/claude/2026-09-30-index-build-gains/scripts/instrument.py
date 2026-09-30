"""Timers and counters for the scaling test (scratch only): CountUniqueKmers split into its all-pairs
flex loop, the rest of the scan and the table output; queries and scanned flex cells of the
uniqueness check. Usage: instrument.py <source dir>"""
import sys

src = sys.argv[1]


def edit(path, pairs):
    path = f"{src}/{path}"
    text = open(path).read()
    for old, new in pairs:
        if text.count(old) != 1:
            sys.exit(f"{path}: {old[:60]!r} found {text.count(old)} times")
        text = text.replace(old, new)
    open(path, "w").write(text)


edit("src/Hash/Seedmap.h", [
    ("#include <atomic>", "#include <atomic>\n#include <chrono>"),
    ("            std::vector<size_t> closest_flex;\n",
     "            std::vector<size_t> closest_flex;\n"
     "            auto I_t0 = std::chrono::steady_clock::now(); std::chrono::nanoseconds I_flex{0};\n"
     "            uint64_t I_pairs = 0, I_flex_keys = 0, I_max_n = 0, I_keys = 0;\n"),
    ("                        auto idx = 0;\n                        auto max_sim = 0;\n",
     "                        auto idx = 0;\n                        auto max_sim = 0;\n"
     "                        auto I_tf = std::chrono::steady_clock::now(); I_pairs += uint64_t(entries) * entries;\n"
     "                        I_flex_keys++; I_max_n = std::max<uint64_t>(I_max_n, entries);\n"),
    ("                            closest_flex[idx++] = m_flex_k - max;\n                        }\n",
     "                            closest_flex[idx++] = m_flex_k - max;\n                        }\n"
     "                        I_flex += std::chrono::steady_clock::now() - I_tf;\n"),
    ("                    bool has_flex_block = key_value_block_size >= m_flex_threshold;\n",
     "                    bool has_flex_block = key_value_block_size >= m_flex_threshold;\n                    I_keys++;\n"),
    ("            for (auto [key, short_uniques] : short_unique_kmers) {",
     "            auto I_t1 = std::chrono::steady_clock::now();\n"
     "            for (auto [key, short_uniques] : short_unique_kmers) {"),
    ('            std::cout << "ShortUniques:    " << short_uniques << std::endl;',
     '            auto I_t2 = std::chrono::steady_clock::now();\n'
     '            std::cerr << "STAT unique_kmers scan_s=" << std::chrono::duration<double>(I_t1 - I_t0).count()\n'
     '                      << " flex_pairs_s=" << std::chrono::duration<double>(I_flex).count()\n'
     '                      << " table_s=" << std::chrono::duration<double>(I_t2 - I_t1).count()\n'
     '                      << " values=" << (short_uniques + long_uniques + non_uniques) << " keys=" << I_keys\n'
     '                      << " flex_keys=" << I_flex_keys << " pairs=" << I_pairs << " max_n=" << I_max_n\n'
     '                      << " genes=" << short_unique_kmers.size() << std::endl;\n'
     '            std::cout << "ShortUniques:    " << short_uniques << std::endl;'),
])

edit("src/Hash/KmerLookup.h", [
    ("        std::vector<uint16_t> flex_vector;\n",
     "        std::vector<uint16_t> flex_vector;\n    public:\n        size_t I_queries = 0, I_scanned = 0, I_single = 0;\n    private:\n"),
    ("        inline void GetFlex(size_t &kmer, std::vector<ValueEntry*>& max_sim_entries, uint32_t& max_similarity) {\n"
     "            m_sm.Get(kmer, m_entry_begin, m_entry_end, m_flex_begin, m_flex_end);\n",
     "        inline void GetFlex(size_t &kmer, std::vector<ValueEntry*>& max_sim_entries, uint32_t& max_similarity) {\n"
     "            m_sm.Get(kmer, m_entry_begin, m_entry_end, m_flex_begin, m_flex_end);\n"
     "            I_queries++; if (m_entry_begin != nullptr && m_flex_begin != nullptr) I_scanned += m_flex_end - m_flex_begin;\n"
     "            else if (m_entry_begin != nullptr && m_entry_end - m_entry_begin == 1) I_single++;\n"),
])

edit("src/Build.h", [
    ("                    single->SetFlagNonUnique();  // atomic\n",
     "                    single->SetFlagNonUnique();  // atomic\n"),
    ("    #pragma omp critical(statistics)\n        statistics.Join(thread_statistics);\n    }\n\n"
     "        std::cout << \"Save unique kmer info",
     "    #pragma omp critical(statistics)\n        {\n        statistics.Join(thread_statistics);\n"
     "        std::cout << \"STAT uniqueness thread queries=\" << lookup.I_queries << \" scanned=\" << lookup.I_scanned\n"
     "                  << \" single=\" << lookup.I_single << \" records=\" << thread_statistics.reads << std::endl;\n"
     "        }\n    }\n\n"
     "        std::cout << \"Save unique kmer info"),
])
print("instrumented", src)
