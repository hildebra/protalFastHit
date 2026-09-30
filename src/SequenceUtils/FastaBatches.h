#pragma once

// A FASTA file in batches of whole records, for the index build's parallel passes (Build.h): one
// thread reads raw bytes and cuts them at record starts, and the thread that takes a batch parses
// its records (ForEachFastaRecord), with what SeqReader gives for FASTA.

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <istream>
#include <string>
#include <string_view>

namespace protal {

    class FastaBatches {
    public:
        explicit FastaBatches(std::istream& is) : m_is(is) {}

        // The next whole records into batch: those that start within its first `bytes` (so about
        // `bytes`, more by the rest of the last record); false at the end of the input or on an
        // error (Error()).
        bool Next(std::string& batch, size_t bytes) {
            batch.swap(m_rest);
            m_rest.clear();
            if (!m_error.empty()) return false;
            bytes = std::max<size_t>(bytes, 1);
            size_t want = bytes;
            while (true) {
                if (!m_end && batch.size() < want) Read(batch, want - batch.size());
                if (batch.empty()) return false;
                if (!m_started) {
                    m_started = true;
                    if (batch[0] != '>') {
                        m_error = "a FASTA file starts with '>', this one with '" + std::string(1, batch[0]) + "'";
                        return false;
                    }
                }
                // Cut before the first record that starts at or after `bytes` (a line end followed by
                // '>'): the batch holds `bytes` and the rest of the record there.
                if (batch.size() > bytes) {
                    size_t const cut = batch.find("\n>", bytes - 1);
                    if (cut != std::string::npos) {
                        m_rest.assign(batch, cut + 1, std::string::npos);
                        batch.resize(cut + 1);
                        return true;
                    }
                }
                if (m_end) return true;  // the rest of the file, whole records
                want = std::max(want, batch.size()) * 2;  // that record goes on: read on
            }
        }

        std::string const& Error() const { return m_error; }

    private:
        void Read(std::string& batch, size_t size) {
            size_t const old = batch.size();
            batch.resize(old + size);
            m_is.read(batch.data() + old, static_cast<std::streamsize>(size));
            size_t const got = static_cast<size_t>(m_is.gcount());
            batch.resize(old + got);
            if (got < size) m_end = true;
        }

        std::istream& m_is;
        std::string m_rest;  // read after the last whole record of the previous batch
        bool m_end = false;
        bool m_started = false;
        std::string m_error;
    };

    // The line without its trailing white space (as BufferedFastxReader's StripString).
    inline std::string_view StrippedLine(std::string_view line) {
        while (!line.empty() && std::isspace(static_cast<unsigned char>(line.back()))) line.remove_suffix(1);
        return line;
    }

    // Calls f(header, sequence) for each record of batch (whole records, as FastaBatches gives them):
    // the header line, and the sequence's lines joined, each without trailing white space, as
    // BufferedFastxReader reads FASTA. A one-line sequence is a view into batch; several lines are
    // joined in scratch.
    template<typename F>
    void ForEachFastaRecord(std::string_view batch, std::string& scratch, F&& f) {
        size_t pos = 0;
        size_t const n = batch.size();
        auto line_end = [&](size_t from) {
            size_t const end = batch.find('\n', from);
            return end == std::string_view::npos ? n : end;
        };
        while (pos < n) {
            size_t end = line_end(pos);
            std::string_view const header = StrippedLine(batch.substr(pos, end - pos));
            pos = end < n ? end + 1 : n;
            std::string_view single;
            size_t lines = 0;
            while (pos < n && batch[pos] != '>') {
                end = line_end(pos);
                std::string_view const line = StrippedLine(batch.substr(pos, end - pos));
                if (lines == 0) {
                    single = line;
                } else {
                    if (lines == 1) scratch.assign(single);
                    scratch.append(line);
                }
                lines++;
                pos = end < n ? end + 1 : n;
            }
            f(header, lines <= 1 ? single : std::string_view(scratch));
        }
    }
}
