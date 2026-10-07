// SPDX-License-Identifier: GPL-2.0-only
#pragma once

// What simulate_metagenomes's read simulators share (LongReadSimulator: long and Ultima reads; IlluminaSimulator:
// Illumina paired-end reads): a random number generator, genomes and a host genome to draw from, and the pipeline
// that makes a run's samples on a pool of threads.
//
// The pipeline (Run): a Job says what a sample's reads are, in rounds of work items (some reads of one genome; a
// sample of long reads plans another round while its bases fall short). Items run on the threads, each with a random
// stream of its own, and their pieces (compressed FASTQ: one zstd frame or BGZF blocks per output file) are appended
// to the sample's files in the items' order: the files do not depend on the threads. Several samples are made at once
// when a sample leaves threads idle; the earlier samples come first. A genome is read once for all the items of one
// round that need it. Outputs that are named pipes (FIFOs) are streamed: one sample at a time, in the samples' order,
// the outputs opened in their order (R1 before R2) and written in place, the pieces of R1 and R2 one after the other
// (each a few hundred kB of FASTQ at most), so that a reader taking R1 and R2 in step (protal) never waits on one
// while the other is full.

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace protal::sim {

// xoshiro256** seeded through splitmix64: small, fast, and the same everywhere.
class LongRng {
public:
    using result_type = std::uint64_t;
    explicit LongRng(std::uint64_t seed);
    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return ~result_type{0}; }
    result_type operator()();
    double Uniform();                      // [0, 1)
    std::uint64_t Below(std::uint64_t n);  // [0, n), n > 0
    double Normal();                       // N(0, 1)

private:
    std::uint64_t m_s[4];
    bool m_has_spare = false;
    double m_spare = 0.0;
};

// A seed made of two numbers (splitmix64 of their mix), for the random streams of work items.
std::uint64_t MixSeed(std::uint64_t a, std::uint64_t b);

// Reverse complement as the collector's COMPLEMENT: ACGTN complemented, other letters kept.
void ReverseComplement(std::string& seq);

// A genome's contigs of at least `min_length` bases, upper case, their names (the first word of the header) and their
// cumulative start ranges (length - min_length + 1 each).
struct Contigs {
    std::vector<std::string> names, seqs;
    std::vector<std::uint64_t> starts;
    std::uint32_t min_length = 1;
};

// Throws (naming `name`) if the file cannot be read or holds no contig of min_length bases.
std::shared_ptr<Contigs const> LoadContigs(std::string const& name, std::filesystem::path const& fasta,
                                           std::uint32_t min_length);

// A host genome prepared by scenarios.prepare_host: host.seq (every contig one after the other) by memory map, and
// host.json's contigs ([name, offset, length], ...).
class Host {
public:
    explicit Host(std::filesystem::path const& folder);
    ~Host();
    Host(Host const&) = delete;
    Host& operator=(Host const&) = delete;

    // As scenarios.Host.draw: `length` bases from a random place (shorter at a contig's end, but at least
    // min(length, min_length)), on the forward strand; one too short or with more than 10% N drawn again, up to 50
    // times.
    std::string Draw(LongRng& rng, std::uint32_t length, std::uint32_t min_length = 100) const;

private:
    std::vector<std::uint64_t> m_offsets, m_lengths, m_ends;
    std::uint64_t m_bases = 0;
    int m_fd = -1;
    std::size_t m_size = 0;
    char const* m_data = nullptr;
};

namespace pipeline {

    // How an output is written, by its name: .zst (zstd) or .gz (BGZF).
    enum class Packing { Zstd, Bgzf };
    Packing PackingOf(std::filesystem::path const& out);
    // A piece of an output: one zstd frame (level 3, with a checksum) or BGZF blocks (the end-of-file block comes when
    // the file is closed); concatenated, they are the file.
    std::string Pack(std::string const& data, Packing packing);

    struct Item {
        int genome = -1;                            // the sample's genome it reads (Job::Fasta), or -1
        std::uint32_t round = 0, part = 0, parts = 1;  // parts: the items of its genome in its round
        std::uint64_t first_read = 0, reads = 0;    // the sample's reads before this item's, and its reads
        std::vector<std::uint32_t> lengths;         // long reads: each read's template length
    };

    struct Piece {
        std::string bytes[2];  // packed, per output file
        std::uint64_t reads = 0, template_bases = 0, read_bases = 0, errors = 0;
    };

    struct Totals {
        std::uint64_t reads = 0, template_bases = 0, read_bases = 0, errors = 0;
        std::uint32_t rounds = 0;  // rounds planned
    };

    class Job {
    public:
        virtual ~Job() = default;
        virtual std::size_t Samples() const = 0;
        // The sample's output files (1 or 2: R1 and R2), .fq.zst or .fq.gz; regular files are written as
        // name.partial and renamed once whole, named pipes in place.
        virtual std::vector<std::filesystem::path> Outputs(std::size_t sample) const = 0;
        // The items of the sample's next round (their round: so_far.rounds), given what its earlier rounds made; none
        // once it is done. Called under the pipeline's lock.
        virtual std::vector<Item> Plan(std::size_t sample, Totals const& so_far) = 0;
        // The FASTA of the sample's genome `genome`, or empty if its items need none (a host).
        virtual std::filesystem::path Fasta(std::size_t sample, int genome) const = 0;
        virtual std::string GenomeName(std::size_t sample, int genome) const = 0;
        virtual std::uint32_t MinContigLength() const = 0;
        // An item's reads, packed (Pack) per output; contigs: its genome's, or null. Called without the lock.
        virtual Piece Make(std::size_t sample, Item const& item, Contigs const* contigs) const = 0;
    };

    // The job's samples on `threads` threads; their totals in the samples' order. Throws on any failure (a genome that
    // cannot be read, a file that cannot be written); the samples written by then stay, the others leave no file.
    std::vector<Totals> Run(Job& job, int threads);

}  // namespace pipeline
}  // namespace protal::sim
