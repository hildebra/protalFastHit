// SPDX-License-Identifier: GPL-2.0-only
#pragma once

// The reference genomes of simulate_metagenomes's read simulators (ReadPipeline's Contigs), read in bulk:
//
// - ReadWholeFile: a file's bytes, inflated whole if gzip (any members, BGZF too) or zstd, told by its first bytes as
//   ThreadedGzStream tells them; no thread and no block buffers per file.
// - ParseFasta: its records as the simulators always took them (a header's first word, the sequence lines without white
//   space, upper case), line by line with memchr and whole lines copied rather than base by base.
// - A genome store (--genome_store DIR): each FASTA decoded once into DIR/<hash of its path>.g2b (GenomeFile), its
//   bases 2-bit and every other character as runs, and memory-mapped when read again, by another sample, another
//   round of long reads or another simulate_metagenomes run (the page cache shares it between processes). A file of
//   the store names the FASTA it was made from, with its size and modification time: a FASTA that changed is decoded
//   again. Reads made from the store are those made from the FASTA, byte for byte.

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace protal::sim {

// Throws (naming `path`) if the file cannot be read, or is truncated or corrupt.
std::string ReadWholeFile(std::filesystem::path const& path);

struct FastaRecords {
    std::vector<std::string> names, seqs;
};

// A record per line starting with '>' (its name: the header up to a space, a tab or a carriage return); its sequence:
// the following lines' characters other than white space (isspace in the C locale), a-z upper-cased; lines before the
// first header are ignored.
FastaRecords ParseFasta(std::string_view text);

// A genome of the store: contigs (names, lengths), bases A, C, G, T in 2 bits each, other characters as runs.
class GenomeFile {
public:
    // The store file of `fasta` in `store`, mapped, if it is there and was made from the FASTA as it is now; else null.
    static std::shared_ptr<GenomeFile const> Open(std::filesystem::path const& store, std::filesystem::path const& fasta);
    // Writes the store file of `fasta` (made from `records`) as a temporary file renamed into place, so that concurrent
    // writers and readers never see half a file. False (and `error`) if it cannot be written.
    static bool Write(std::filesystem::path const& store, std::filesystem::path const& fasta, FastaRecords const& records,
                      std::string& error);
    // Where the store keeps `fasta`'s file.
    static std::filesystem::path PathOf(std::filesystem::path const& store, std::filesystem::path const& fasta);

    ~GenomeFile();
    GenomeFile(GenomeFile const&) = delete;
    GenomeFile& operator=(GenomeFile const&) = delete;

    std::size_t Contigs() const { return m_contigs; }
    std::string_view Name(std::size_t contig) const;
    std::uint64_t Length(std::size_t contig) const;
    // As out.assign(sequence, start, length) of the contig's sequence (start at most its length).
    void Extract(std::size_t contig, std::uint64_t start, std::uint64_t length, std::string& out) const;

private:
    GenomeFile() = default;
    struct Entry;
    struct Run;
    Entry const& EntryOf(std::size_t contig) const;

    char const* m_data = nullptr;
    std::size_t m_size = 0;
    std::size_t m_contigs = 0;
    Entry const* m_entries = nullptr;
    char const* m_names = nullptr;
    Run const* m_runs = nullptr;
};

}  // namespace protal::sim
