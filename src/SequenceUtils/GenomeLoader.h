//
// Created by fritsche on 15/07/22.
//

#pragma once

#include <Constants.h>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <sparse_map.h>
#include <string>
#include <fstream>
#include <err.h>
#include <KmerUtils.h>


#include <sparse_set.h>

#include "Utilities.h"
#include "Zstd.h"
#include <sysexits.h>

#include "Benchmark.h"

namespace protal {
    class Gene {
        static const size_t DEFAULT = UINT64_MAX;
    private:
        size_t m_id = 0;
        std::string m_sequence = "";
        size_t m_start_byte;
        size_t m_length = DEFAULT;
        std::ifstream* m_is = nullptr;  // nullptr: compressed reference, genes only come from preloading

        size_t m_short_unique = 0;
        size_t m_long_unique = 0;
        size_t m_long_super_unique = 0;
        size_t m_total_kmers = 0;

        // Reads the gene's bytes from reference.fna at the offsets reference.map gives. A failed read
        // stops protal: the stream is shared by all genes, and leaving it failed would silently turn
        // every gene read after it into NUL bytes. Sequences are uppercased (the k-mer and alignment
        // code only knows A, C, G, T).
        void Load(std::string& into, size_t start_byte, size_t length) {
            if(m_is && m_is->is_open())
            {
                m_is->seekg(start_byte);
                into.resize(length);
                m_is->read(&into[0], length);
                if (!*m_is || static_cast<size_t>(m_is->gcount()) != length) {
                    std::cerr << "Cannot read gene " << m_id << " (bytes " << start_byte << "-" << start_byte + length
                              << ") from reference.fna: reference.map does not match the file" << std::endl;
                    exit(8);
                }
                Uppercase(into);
            } else if (!m_is && length > 0 && length != DEFAULT) {
                errx(EX_SOFTWARE, "Gene %zu cannot be loaded on its own from a compressed reference (reference.fna.zst); "
                                  "the reference must be preloaded.", m_id);
            }
        };

    public:
        // The k-mer and alignment code only knows A, C, G, T: lowercase bases would encode as A.
        static void Uppercase(std::string& sequence) {
            std::transform(sequence.begin(), sequence.end(), sequence.begin(), [](unsigned char c) { return std::toupper(c); });
        }

        Gene(){};

        void Set(size_t id, size_t start_byte, size_t length, std::ifstream* is) {
            m_id = id;
            m_start_byte = start_byte;
            m_length = length;
            m_is = is;
        }

        size_t GetStartByte() const {
            return m_start_byte;
        }

        void SetSequence(std::string&& sequence) {
            m_sequence = std::move(sequence);
        }

        bool IsSet() const {
            return (m_id != 0 || m_length > 0) && m_length != DEFAULT;
        }

        bool IsLoaded() const {
            return !m_sequence.empty();
        }

        void Load() {
            if (!IsLoaded()) {
                Load(m_sequence, m_start_byte, m_length);
            }
        };

        void LoadOMP() {
#pragma omp critical(genome_loader)
            if (!IsLoaded()) {
                Load(m_sequence, m_start_byte, m_length);
            }
        };

        size_t GetLength() const {
            return m_length;
        }

        const std::string& Sequence() const {
            return m_sequence;
        }
        const size_t GetId() const {
            return m_id;
        }

        void SetUniqueValues(size_t short_unique, size_t long_unique, size_t long_super_unique, size_t total_kmers) {
            m_short_unique = short_unique;
            m_long_unique = long_unique;
            m_long_super_unique = long_super_unique;
            m_total_kmers = total_kmers;
        }

        [[nodiscard]] std::tuple<size_t, size_t, size_t, size_t> GetUniqueKmerCounts() const {
            return { m_short_unique, m_long_unique, m_long_super_unique, m_total_kmers };
        }

        bool HasShortUniques(const size_t threshold=0) const {
            return m_short_unique > threshold;
        }
        bool HasLongUniques(const size_t threshold=0) const {
            return m_long_unique > threshold;
        }
        bool HasLongSuperUniques(const size_t threshold=0) const {
            return m_long_super_unique > threshold;
        }

        double UniqueRate() const {
            return m_long_unique/static_cast<double>(m_total_kmers);
        }

        double SuperUniqueRate() const {
            return m_long_super_unique/static_cast<double>(m_total_kmers);
        }
    };


    class Genome {
    public:
        using GeneKey = size_t;
        using GenomeKey = size_t;
        using GeneList = std::vector<Gene>;

    private:
        using GeneID = uint32_t;
        GeneList m_genes;
        GenomeKey m_key;
        tsl::sparse_set<GeneID> m_hittable_genes;
        bool m_is_loaded = false;

        size_t m_short_unique = 0;
        size_t m_long_unique = 0;
        size_t m_long_super_unique = 0;
        size_t m_total_kmers = 0;



        const size_t GeneKeyToIndex(GeneKey const& key) const {
            return key - 1;
        }
    public:

//        Genome() {};
        explicit Genome(GenomeKey key) : m_key(key) {};

        GenomeKey GetKey() {
            return m_key;
        }

        void SetUniqueValues() {
            for (auto& gene : m_genes) {
                auto [su, lu, lsu, total] = gene.GetUniqueKmerCounts();
                m_short_unique += su;
                m_long_unique += lu;
                m_long_super_unique += lsu;
                m_total_kmers += total;
            }
        }

        [[nodiscard]] std::tuple<size_t, size_t, size_t, size_t> GetUniqueKmerCounts() const {
            return { m_short_unique, m_long_unique, m_long_super_unique, m_total_kmers };
        }

        void AddGene(GeneKey key, size_t id, size_t start_byte, size_t length, std::ifstream* is) {
            auto index = GeneKeyToIndex(key);
            if (index >= m_genes.size()) {
                m_genes.resize(index+1);
            }
            m_genes[index].Set(key, start_byte, length, is);
        }

        void AddHittableGene(GeneID geneid) {
            m_hittable_genes.insert(geneid);
        }

        bool IsGeneHittable(GeneID geneid) const {
            return m_hittable_genes.empty() ? true : m_hittable_genes.contains(geneid);
        }

        size_t GenesWithShortUniques(size_t threshold = 0) {
            return std::count_if(m_genes.begin(), m_genes.end(), [threshold](Gene const& gene) {
                return gene.HasShortUniques(threshold);
            });
        }

        size_t GenesWithLongUniques(size_t threshold = 0) const {
            return std::count_if(m_genes.begin(), m_genes.end(), [threshold](Gene const& gene) {
                return gene.HasLongUniques(threshold);
            });
        }

        size_t GenesWithLongSuperUniques(size_t threshold = 0) const {
            return std::count_if(m_genes.begin(), m_genes.end(), [threshold](Gene const& gene) {
                return gene.HasLongSuperUniques(threshold);
            });
        }

        std::vector<uint32_t> GetHittableGenes() {
            std::vector<uint32_t> genes;
            if (m_hittable_genes.empty()) {
                genes.reserve(m_genes.size());
                for (auto const& gene : m_genes) {
                    if (gene.IsSet()) {
                        genes.emplace_back(gene.GetId());
                    }
                }
            } else {
                genes.reserve(m_hittable_genes.size());
                for (auto const& gene_id : m_hittable_genes) {
                    genes.emplace_back(gene_id);
                }
                std::sort(genes.begin(), genes.end());
            }
            return genes;
        }

        void LoadGene(GeneKey key) {
            m_genes[GeneKeyToIndex(key)].Load();
        };

        bool ValidGene(GeneKey key) {
            return GeneKeyToIndex(key) < m_genes.size();
        }

        bool HasGene(GeneKey key) const {
            return key > 0 && GeneKeyToIndex(key) < m_genes.size() && m_genes[GeneKeyToIndex(key)].IsSet();
        }

        Gene& GetGene(GeneKey key) {
            return m_genes.at(GeneKeyToIndex(key));
        }

        const GeneList GetGeneList() const {
            return m_genes;
        }

        GeneList& Genes() {
            return m_genes;
        }

        void MarkLoaded() {
            m_is_loaded = true;
        }

        size_t GeneNum() const {
            return m_hittable_genes.empty() ? std::count_if(m_genes.begin(), m_genes.end(), [](Gene const& gene) {
                return gene.IsSet();
            }) : m_hittable_genes.size();
        }

        void LoadGenome() {
            std::for_each(m_genes.begin(), m_genes.end(), [](Gene& gene){ gene.Load(); });
            m_is_loaded = true;
        };

        bool IsLoaded() const {
            return m_is_loaded;
//            return std::any_of(m_genes.begin(), m_genes.end(), [](Gene const& gene){ return gene.IsLoaded(); });
        }

        void LoadGenomeOMP() {
#pragma omp critical(genome_loader)
            {
                if (!IsLoaded()) {
                    std::for_each(m_genes.begin(), m_genes.end(), [](Gene &gene) { gene.Load(); });
                    m_is_loaded = true;
                }
            }
        };

        Gene& GetGeneOMP(GeneKey key) {
            if (!m_is_loaded) {
#pragma omp critical(genome_loader)
                if (!m_is_loaded) {
                    std::for_each(m_genes.begin(), m_genes.end(), [](Gene &gene) { gene.Load(); });
                    m_is_loaded = true;
                }
            }
            return m_genes.at(GeneKeyToIndex(key));
        };
    };


    class GenomeLoader {
        using GenomeKey = Genome::GenomeKey;
        using GenomeMap = tsl::sparse_map<GenomeKey, Genome>;
        using GeneKey = Genome::GeneKey;


        std::string m_path;
        std::string m_genome_map;
        bool m_compressed = false;  // reference.fna.zst: genes are only read by LoadAllGenomes
        std::ifstream m_is;
        GenomeMap m_genomes;

        Genome& AddOrGetGenome(GenomeKey const& key) {
            if (!m_genomes.contains(key)) {
                m_genomes.insert( { key, Genome(key) } );
            }
            return m_genomes.at(key);
        }

        void Open() {
            m_compressed = zstd::IsCompressed(m_path);
            if (!m_compressed) m_is.open(m_path, std::ios::in);
        }

    public:
        // genome_path: reference.fna or a zstd-compressed reference.fna.zst. The byte offsets in
        // genome_map (reference.map) always refer to the uncompressed reference.
        GenomeLoader(std::string genome_path, std::string genome_map) :
                m_path(genome_path),
                m_genome_map(genome_map) {
            Open();
            LoadPositionMap(genome_map);
        };

        GenomeLoader(const GenomeLoader& other) :
                m_path(other.m_path),
                m_genome_map(other.m_genome_map) {
            Open();
            LoadPositionMap(other.m_genome_map);
        }

        bool IsCompressed() const {
            return m_compressed;
        }

        ~GenomeLoader() {
            m_is.close();
        }

        void PrintHittableGenes() {
            for (auto& [gid, _] : m_genomes) {
                auto& genome = m_genomes.at(gid);
                auto hg = genome.GetHittableGenes();
                std::string hgstr = "";
                for (auto gene : hg) {
                    hgstr += std::to_string(gene) + ',';
                }
                std::cout << "Hittable\t" << gid << '\t' << hg.size() << '\t' << hgstr << std::endl;
            }
        }

        // unique_kmers.tsv (written by --build): taxid, gene id, then counts and rates of short,
        // long and long-super unique k-mers, and the gene's k-mer total. Every line must name a
        // gene of reference.map.
        void LoadUniqueKmers(std::string const& file) {
            std::ifstream is(file, std::ios::in);
            if (!is) InvalidUniqueKmers(file, 0, "cannot open the file");

            std::vector<std::string> tokens;
            std::string line;
            size_t line_no = 0;
            while (std::getline(is, line)) {
                line_no++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty()) continue;
                Utils::split(tokens, line, "\t");
                if (tokens.size() < 9) {
                    InvalidUniqueKmers(file, line_no, "expected 9 tab-separated columns, found " + std::to_string(tokens.size()));
                }
                uint64_t counts[5];
                size_t const count_columns[5] = { 0, 1, 2, 4, 6 };
                for (int i = 0; i < 5; i++) {
                    auto const& t = tokens[count_columns[i]];
                    if (t.empty() || !std::all_of(t.begin(), t.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); })) {
                        InvalidUniqueKmers(file, line_no, "column " + std::to_string(count_columns[i] + 1) + " is not a non-negative integer");
                    }
                    counts[i] = std::stoull(t);
                }
                auto const& total_str = tokens[8];
                if (total_str.empty() || !std::all_of(total_str.begin(), total_str.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); })) {
                    InvalidUniqueKmers(file, line_no, "column 9 is not a non-negative integer");
                }

                auto taxid = counts[0];
                auto geneid = counts[1];
                auto short_unique = counts[2];
                auto long_unique = counts[3];
                auto long_super_unique = counts[4];
                auto total_kmers = std::stoull(total_str);
                if (!HasGene(taxid, geneid)) {
                    InvalidUniqueKmers(file, line_no, "gene " + std::to_string(taxid) + "_" + std::to_string(geneid) +
                                                      " is not in reference.map (rebuild the database with --build)");
                }

                auto& taxon = m_genomes.at(taxid);
                if (short_unique + long_unique > 0) {
                    taxon.AddHittableGene(geneid);
                }
                auto& gene = taxon.GetGene(geneid);
                gene.SetUniqueValues(short_unique, long_unique, long_super_unique, total_kmers);
            }
            is.close();


            for (auto& tid : m_genomes) {
                m_genomes.at(tid.first).SetUniqueValues();
            }

            // PrintHittableGenes();
        }

        void LoadHittableGenes(std::string const& file) {
            std::ifstream is(file, std::ios::in);

            std::vector<std::string> tokens;
            std::string line;
            while (std::getline(is, line)) {
                Utils::split(tokens, line, "\t");
                auto taxid = std::stoull(tokens[0]);
                auto gene_str = tokens[2];

                auto& taxon = m_genomes.at(taxid);

                Utils::split(tokens, gene_str, ",");
                for (auto const& g : tokens) {
                    auto gid = std::stoul(g);
                    taxon.AddHittableGene(gid);
                }
            }
            is.close();
        }

        size_t GetLoadedGenomeCount() const {
            return std::count_if(m_genomes.begin(), m_genomes.end(), [](auto const& pair) { return pair.second.IsLoaded(); });
        }

        const Genome& GetGenome(GenomeKey const& key) const {
            return m_genomes.at(key);
        }

        bool HasGene(GenomeKey taxid, GeneKey gene) const {
            auto it = m_genomes.find(taxid);
            return it != m_genomes.end() && it->second.HasGene(gene);
        }

        // Length of a gene as reference.map gives it (0 if the gene is not in the map).
        size_t GeneLength(GenomeKey taxid, GeneKey gene) {
            return HasGene(taxid, gene) ? m_genomes.at(taxid).GetGene(gene).GetLength() : 0;
        }

        size_t GeneCount() const {
            size_t n = 0;
            for (auto const& [key, genome] : m_genomes) {
                for (auto const& gene : genome.GetGeneList()) n += gene.IsSet();
            }
            return n;
        }

        GenomeMap& GetGenomeMap() {
            return m_genomes;
        }

        void WriteSamHeader(std::ostream& os=std::cout) {
            os << "@HD\tVN:1.6\n";
            for (auto& [key, genome] : m_genomes) {
                auto& genes = genome.GetGeneList();
                for (auto i = 0; i < genes.size(); i++) {
                    if (genes[i].IsSet()) {

                        os << "@SQ\tSN:" << key << '_' << genes[i].GetId() << '\t' << "LN:" << genes[i].GetLength() << '\n';
                    }
                }
            }
        }

        // Reads every gene in one sequential pass over the reference (raw or compressed), in
        // file order, skipping the header lines between them: large sequential reads instead of
        // one seek per gene, which matters on network storage.
        void LoadAllGenomes() {
            std::vector<GenomeKey> keys;
            for (auto& pair : m_genomes) {
                keys.emplace_back(pair.first);
            }
            std::sort(keys.begin(), keys.end());

            std::vector<Gene*> genes;
            for (auto& key : keys) {
                auto& genome = m_genomes.at(key);
                if (genome.IsLoaded()) continue;
                for (auto& gene : genome.Genes()) {
                    if (gene.IsSet() && !gene.IsLoaded() && gene.GetLength() > 0) genes.emplace_back(&gene);
                }
            }
            std::sort(genes.begin(), genes.end(), [](Gene const* a, Gene const* b) {
                return a->GetStartByte() < b->GetStartByte();
            });

            zstd::InputFile input(m_path);
            if (!input.IsOpen()) {
                std::cerr << "Cannot open the reference " << m_path << std::endl;
                exit(8);
            }
            std::istream& is = input.Stream();
            size_t position = 0;
            for (Gene* gene : genes) {
                if (gene->GetStartByte() < position) {
                    std::cerr << "Invalid reference map " << m_genome_map << ": gene " << gene->GetId() << " at byte "
                              << gene->GetStartByte() << " overlaps the previous gene" << std::endl;
                    exit(8);
                }
                is.ignore(static_cast<std::streamsize>(gene->GetStartByte() - position));
                std::string sequence(gene->GetLength(), '\0');
                is.read(sequence.data(), static_cast<std::streamsize>(sequence.size()));
                if (!is) {
                    std::cerr << "Cannot read gene " << gene->GetId() << " (bytes " << gene->GetStartByte() << "-"
                              << gene->GetStartByte() + gene->GetLength() << ") from " << m_path
                              << ": reference.map does not match the file" << std::endl;
                    exit(8);
                }
                Gene::Uppercase(sequence);
                gene->SetSequence(std::move(sequence));
                position = gene->GetStartByte() + gene->GetLength();
            }
            for (auto& key : keys) {
                m_genomes.at(key).MarkLoaded();
            }
        }

        bool AllGenomesLoaded() {
            for (auto& [id, genome] : m_genomes) {
                if (!genome.IsLoaded()){
                    return false;
                }
            }
            return true;
        }

        Genome& GetGenome(GenomeKey const& key) {
            assert(m_genomes.contains(key));
            if (!m_genomes.contains(key)) {
                std::cout << "Genomes Key: " << key << std::endl;
                exit(10);
            }
            return m_genomes.at(key);
        }


        // reference.map: taxid, gene id, start byte, end byte of the gene's sequence in reference.fna.
        // Every line is checked, so a bad map stops protal here instead of corrupting genes silently.
        void LoadPositionMap(std::string file_path) {
            std::ifstream is(file_path, std::ios::in);
            if (!is) InvalidMap(file_path, 0, "cannot open the file");
            // Offsets refer to the uncompressed reference, also for reference.fna.zst.
            auto const size = zstd::UncompressedSize(m_path);
            if (!size) InvalidMap(file_path, 0, "cannot read the size of " + m_path);
            uint64_t const fna_size = *size;
            constexpr uint64_t max_id = (uint64_t{1} << SEEDMAP_TAXID_BITS) - 1;
            constexpr uint64_t max_gene = (uint64_t{1} << SEEDMAP_GENEID_BITS) - 1;
            constexpr uint64_t max_length = (uint64_t{1} << SEEDMAP_GENE_POS_BITS) - 1;

            std::string line;
            size_t line_no = 0;
            while (std::getline(is, line)) {
                line_no++;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (line.empty()) continue;
                auto tokens = Utils::split(line, "\t");

                if (tokens.size() != 4) {
                    InvalidMap(file_path, line_no, "expected 4 tab-separated columns, found " + std::to_string(tokens.size()));
                }
                uint64_t values[4];
                for (int i = 0; i < 4; i++) {
                    auto const& t = tokens[i];
                    if (t.empty() || !std::all_of(t.begin(), t.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); })) {
                        InvalidMap(file_path, line_no, "column " + std::to_string(i + 1) + " is not a non-negative integer");
                    }
                    values[i] = std::stoull(t);
                }
                GenomeKey genome_id = values[0];
                GeneKey gene_key = values[1];
                size_t start = values[2];
                size_t end = values[3];

                // taxid 0 marks empty index entries and gene ids are 1-based; both are packed into 20 bits.
                if (genome_id == 0 || genome_id > max_id) InvalidMap(file_path, line_no, "taxid must be between 1 and " + std::to_string(max_id));
                if (gene_key == 0 || gene_key > max_gene) InvalidMap(file_path, line_no, "gene id must be between 1 and " + std::to_string(max_gene));
                if (end <= start) InvalidMap(file_path, line_no, "end byte must be after start byte");
                if (end > fna_size) InvalidMap(file_path, line_no, "end byte " + std::to_string(end) + " is past the end of " + m_path + " (" + std::to_string(fna_size) + " bytes)");
                if (end - start > max_length) InvalidMap(file_path, line_no, "gene is longer than " + std::to_string(max_length) + " bases");

                auto& genome = AddOrGetGenome(genome_id);
                if (genome.HasGene(gene_key)) InvalidMap(file_path, line_no, "gene " + std::to_string(genome_id) + "_" + std::to_string(gene_key) + " is listed twice");

                genome.AddGene(gene_key, gene_key, start, end - start, m_compressed ? nullptr : &m_is);
            }
            is.close();
        }

        [[noreturn]] static void InvalidUniqueKmers(std::string const& path, size_t line_no, std::string const& reason) {
            std::cerr << "Invalid unique k-mer file " << path;
            if (line_no > 0) std::cerr << ", line " << line_no;
            std::cerr << ": " << reason << std::endl;
            exit(8);
        }

        [[noreturn]] static void InvalidMap(std::string const& path, size_t line_no, std::string const& reason) {
            std::cerr << "Invalid reference map " << path;
            if (line_no > 0) std::cerr << ", line " << line_no;
            std::cerr << ": " << reason << std::endl;
            exit(8);
        }
    };
}
