// GeneTableFile.h - the binary gene table, gene_table.bin: reference.map and unique_kmers.tsv as GenomeLoader holds
// them, by genome, in a member of the single-file database beside the two text tables. A run from a single-file
// database loads it instead of them: no parsing, every genome's genes added by one thread, and the reference
// fingerprint the index is checked against without reading reference.map again. At GTDB r226 size (14.5M genes) the
// two text tables took 3.4 s of every run on 32 threads, parsing most of it
// (docs/claude/2026-10-04-performance-gtdb-scale).
//
// The table is derived: the packer writes it from the text tables it packs (protal --build, --compress_db, which also
// adds it to a single-file database without one), checked as a run checks them (GenomeLoader); --unpack_db leaves it
// out, and a database folder uses the text tables. A table whose recorded member sizes or reference size differ from
// the database's is not used.
//
// Layout (little-endian, as protal's index):
//   Header   magic "PRGENES1", version, the sizes of the reference.map and unique_kmers.tsv members it was made from
//            (unique_kmers.tsv: 0 if there was none), whether it holds the unique k-mer counts, the ReferenceFingerprint
//            (xxhash64 of reference.map, size of reference.fna), the number of genomes and of genes
//   Genomes  in the order reference.map lists them (first appearance; the loader makes the genomes in that order, as the
//            text loader does): taxid, gene slots (the largest gene id), the index of its first gene, its number of genes
//   Genes    by genome, by gene id: start byte in reference.fna, gene id, length, the four unique k-mer counts
#pragma once

#include <cstdint>
#include <string>
#include <type_traits>

namespace protal::gene_table_file {
    inline const std::string kFileName = "gene_table.bin";
    inline constexpr char kMagic[8] = {'P', 'R', 'G', 'E', 'N', 'E', 'S', '1'};
    inline constexpr uint64_t kVersion = 1;

    struct Header {
        char magic[8];
        uint64_t version;
        uint64_t map_size;      // the reference.map member's (uncompressed) size
        uint64_t unique_size;   // the unique_kmers.tsv member's, 0 if none
        uint64_t has_unique;    // 1: the genes carry their unique k-mer counts (unique_kmers.tsv)
        uint64_t map_hash;      // ReferenceFingerprint::map_hash of reference.map
        uint64_t fna_size;      // ReferenceFingerprint::fna_size: reference.fna's uncompressed size
        uint64_t genomes;
        uint64_t genes;
    };

    struct Genome {
        uint64_t taxid;
        uint64_t slots;  // the largest gene id: the genome's gene list is this long
        uint64_t first;  // its first gene among the genes
        uint64_t count;
    };

    struct Gene {
        uint64_t start;  // byte in reference.fna (uncompressed)
        uint32_t id;
        uint32_t length;
        uint32_t short_unique, long_unique, long_super_unique, total_kmers;
    };

    static_assert(sizeof(Header) == 72 && sizeof(Genome) == 32 && sizeof(Gene) == 32, "gene_table.bin record sizes");
    static_assert(std::is_trivially_copyable_v<Header> && std::is_trivially_copyable_v<Genome> && std::is_trivially_copyable_v<Gene>);

    inline uint64_t ExpectedSize(Header const& h) {
        return sizeof(Header) + h.genomes * sizeof(Genome) + h.genes * sizeof(Gene);
    }
}
