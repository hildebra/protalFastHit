// SPDX-License-Identifier: GPL-2.0-only
#include "ReadPipeline.h"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <exception>
#include <fstream>
#include <future>
#include <iterator>
#include <map>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>
#include <utility>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../IO/Bgzf.h"
#include "../Utilities/Zstd.h"
#include "ThreadedGzStream.h"

namespace fs = std::filesystem;

namespace protal::sim {

// ---- random numbers ---------------------------------------------------------------------------------------------

static std::uint64_t SplitMix64(std::uint64_t& x) {
    std::uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static std::uint64_t Rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

LongRng::LongRng(std::uint64_t seed) {
    for (auto& s : m_s) s = SplitMix64(seed);
}

LongRng::result_type LongRng::operator()() {
    std::uint64_t const result = Rotl(m_s[1] * 5, 7) * 9;
    std::uint64_t const t = m_s[1] << 17;
    m_s[2] ^= m_s[0];
    m_s[3] ^= m_s[1];
    m_s[1] ^= m_s[2];
    m_s[0] ^= m_s[3];
    m_s[2] ^= t;
    m_s[3] = Rotl(m_s[3], 45);
    return result;
}

double LongRng::Uniform() { return static_cast<double>((*this)() >> 11) * 0x1.0p-53; }

std::uint64_t LongRng::Below(std::uint64_t n) {  // Lemire's multiply-and-reject: unbiased
    __uint128_t m = static_cast<__uint128_t>((*this)()) * n;
    auto low = static_cast<std::uint64_t>(m);
    if (low < n) {
        std::uint64_t const threshold = -n % n;
        while (low < threshold) {
            m = static_cast<__uint128_t>((*this)()) * n;
            low = static_cast<std::uint64_t>(m);
        }
    }
    return static_cast<std::uint64_t>(m >> 64);
}

double LongRng::Normal() {  // Marsaglia's polar method
    if (m_has_spare) {
        m_has_spare = false;
        return m_spare;
    }
    double u, v, s;
    do {
        u = 2.0 * Uniform() - 1.0;
        v = 2.0 * Uniform() - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    double const f = std::sqrt(-2.0 * std::log(s) / s);
    m_spare = v * f;
    m_has_spare = true;
    return u * f;
}

std::uint64_t MixSeed(std::uint64_t a, std::uint64_t b) {
    std::uint64_t x = a ^ Rotl(b, 32) ^ 0x2545f4914f6cdd1dULL;
    SplitMix64(x);
    return SplitMix64(x);
}

void ReverseComplement(std::string& seq) {
    std::reverse(seq.begin(), seq.end());
    for (char& c : seq) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            default: break;  // N and others as they are
        }
    }
}

// ---- genomes -----------------------------------------------------------------------------------------------------

std::shared_ptr<Contigs const> LoadContigs(std::string const& name, fs::path const& fasta, std::uint32_t min_length) {
    protal::ThreadedGzIstream in(fasta.c_str());
    if (!in.rdbuf()->is_open()) throw std::runtime_error(name + ": cannot read " + fasta.string());
    auto contigs = std::make_shared<Contigs>();
    contigs->min_length = std::max<std::uint32_t>(1, min_length);
    std::string line, seq, contig;
    bool in_record = false;
    auto close = [&] {
        if (in_record && seq.size() >= contigs->min_length) {
            contigs->starts.push_back((contigs->starts.empty() ? 0 : contigs->starts.back()) + seq.size() -
                                      contigs->min_length + 1);
            contigs->seqs.push_back(std::move(seq));
            contigs->names.push_back(contig);
        }
        seq.clear();
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            close();
            in_record = true;
            auto const end = line.find_first_of(" \t\r", 1);
            contig = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
            continue;
        }
        if (!in_record) continue;
        for (char c : line) {
            if (!std::isspace(static_cast<unsigned char>(c))) seq += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        }
    }
    if (in.rdbuf()->read_failed()) {
        throw std::runtime_error(name + ": " + fasta.string() + " is truncated or corrupt (" +
                                 in.rdbuf()->read_error_message() + ")");
    }
    close();
    if (contigs->seqs.empty()) {
        throw std::runtime_error(name + ": no sequence of " + std::to_string(contigs->min_length) + " bases or more in " +
                                 fasta.string());
    }
    return contigs;
}

Host::Host(fs::path const& folder) {
    std::ifstream in(folder / "host.json");
    if (!in) throw std::runtime_error("cannot read the host index " + (folder / "host.json").string());
    std::string const json((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    auto fail = [&](std::string const& why) {
        return std::runtime_error("host index " + (folder / "host.json").string() + ": " + why);
    };
    std::size_t at = json.find("\"contigs\"");
    if (at == std::string::npos || (at = json.find('[', at)) == std::string::npos) throw fail("no contigs");
    ++at;
    auto skip = [&] {
        while (at < json.size() && std::isspace(static_cast<unsigned char>(json[at]))) ++at;
    };
    auto number = [&]() -> std::uint64_t {
        skip();
        std::size_t used = 0;
        std::uint64_t const v = std::stoull(json.substr(at, 24), &used);
        at += used;
        return v;
    };
    auto expect = [&](char c) {
        skip();
        if (at >= json.size() || json[at] != c) throw fail(std::string("expected '") + c + "'");
        ++at;
    };
    for (;;) {
        skip();
        if (at < json.size() && json[at] == ']') break;
        expect('[');
        expect('"');
        while (at < json.size() && json[at] != '"') at += json[at] == '\\' ? 2 : 1;  // the name, not needed
        expect('"');
        expect(',');
        std::uint64_t const offset = number();
        expect(',');
        std::uint64_t const length = number();
        expect(']');
        m_offsets.push_back(offset);
        m_lengths.push_back(length);
        m_bases += length;
        m_ends.push_back(m_bases);
        skip();
        if (at < json.size() && json[at] == ',') ++at;
    }
    if (m_ends.empty()) throw fail("no contigs");
    fs::path const seq = folder / "host.seq";
    m_fd = ::open(seq.c_str(), O_RDONLY);
    if (m_fd < 0) throw std::runtime_error("cannot read " + seq.string() + ": " + std::strerror(errno));
    struct stat st {};
    if (::fstat(m_fd, &st) != 0) throw std::runtime_error("cannot read " + seq.string() + ": " + std::strerror(errno));
    m_size = static_cast<std::size_t>(st.st_size);
    for (std::size_t k = 0; k < m_offsets.size(); ++k) {
        if (m_offsets[k] + m_lengths[k] > m_size) throw std::runtime_error(seq.string() + " is shorter than its index says");
    }
    void* map = ::mmap(nullptr, m_size, PROT_READ, MAP_SHARED, m_fd, 0);
    if (map == MAP_FAILED) throw std::runtime_error("cannot map " + seq.string() + ": " + std::strerror(errno));
    m_data = static_cast<char const*>(map);
}

Host::~Host() {
    if (m_data) ::munmap(const_cast<char*>(m_data), m_size);
    if (m_fd >= 0) ::close(m_fd);
}

std::string Host::Draw(LongRng& rng, std::uint32_t length, std::uint32_t min_length) const {
    std::uint64_t const least = std::min<std::uint64_t>(length, min_length);
    std::string seq;
    for (int attempt = 0; attempt < 50; ++attempt) {
        std::uint64_t const at = rng.Below(m_bases);
        std::size_t const k = std::upper_bound(m_ends.begin(), m_ends.end(), at) - m_ends.begin();
        std::uint64_t const start = at - (m_ends[k] - m_lengths[k]);
        std::uint64_t const end = std::min<std::uint64_t>(m_lengths[k], start + length);
        seq.assign(m_data + m_offsets[k] + start, end - start);
        if (seq.size() >= least && std::count(seq.begin(), seq.end(), 'N') <= 0.1 * static_cast<double>(seq.size())) break;
    }
    return seq;
}

// ---- the pipeline ------------------------------------------------------------------------------------------------

namespace pipeline {

Packing PackingOf(fs::path const& out) {
    std::string const name = out.filename().string();
    auto ends = [&](std::string const& suffix) {
        return name.size() > suffix.size() && name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0;
    };
    if (ends(".zst")) return Packing::Zstd;
    if (ends(".gz")) return Packing::Bgzf;
    throw std::runtime_error(out.string() + ": reads go to .fq.zst (zstd) or .fq.gz (BGZF)");
}

std::string Pack(std::string const& data, Packing packing) {
    std::string out;
    if (packing == Packing::Bgzf) {
        if (!protal::bgzf::Compress(data.data(), data.size(), out)) throw std::runtime_error("BGZF compression failed");
        return out;
    }
    struct Free {
        void operator()(ZSTD_CCtx* c) const { ZSTD_freeCCtx(c); }
    };
    thread_local std::unique_ptr<ZSTD_CCtx, Free> cctx;
    if (!cctx) {
        std::string error;
        cctx.reset(protal::zstd::MakeCCtx(protal::zstd::Params{3, 0, 1, 0}, error));
        if (!cctx) throw std::runtime_error(error);
    }
    out.resize(ZSTD_compressBound(data.size()));
    std::size_t const size = ZSTD_compress2(cctx.get(), out.data(), out.size(), data.data(), data.size());
    if (ZSTD_isError(size)) throw std::runtime_error(std::string("zstd compression failed: ") + ZSTD_getErrorName(size));
    out.resize(size);
    return out;
}

namespace {

bool IsFifo(fs::path const& path) {
    struct stat st {};
    return ::stat(path.c_str(), &st) == 0 && S_ISFIFO(st.st_mode);
}

struct GenomeSlot {
    std::shared_future<std::shared_ptr<Contigs const>> loaded;
    std::uint32_t parts_left = 0;
};

struct Output {
    fs::path path, written;  // written: path.partial, or the pipe itself
    Packing packing = Packing::Zstd;
    bool fifo = false;
    std::FILE* file = nullptr;
};

struct Stream {
    std::size_t index = 0;
    std::vector<Output> outputs;
    std::vector<Item> items;
    std::vector<std::optional<Piece>> pieces;
    std::size_t next_issue = 0, next_write = 0, processed = 0, round_end = 0;
    bool planned_all = false, writing = false, closed = false;
    std::map<std::pair<std::uint32_t, int>, GenomeSlot> genomes;
    Totals totals;
};

class Engine {
public:
    Engine(Job& job, int threads) : m_job(job), m_streams(job.Samples()) {
        for (std::size_t s = 0; s < m_streams.size(); ++s) {
            for (auto const& path : job.Outputs(s)) {
                PackingOf(path);  // fails before any work
                m_streaming = m_streaming || IsFifo(path);
            }
        }
        m_threads = std::max(1, threads);
        m_cap = std::max<std::size_t>(4, 3 * static_cast<std::size_t>(m_threads));
    }

    ~Engine() {
        for (auto& s : m_streams) {
            for (auto& out : s.outputs) {
                if (!out.file) continue;
                std::fclose(out.file);
                if (!out.fifo) {
                    std::error_code ec;
                    fs::remove(out.written, ec);
                }
            }
        }
    }

    std::vector<Totals> Run() {
        std::vector<std::thread> pool;
        for (int t = 1; t < m_threads; ++t) pool.emplace_back([this] { Work(); });
        Work();
        for (auto& thread : pool) thread.join();
        if (m_error) std::rethrow_exception(m_error);
        std::vector<Totals> totals;
        for (auto const& s : m_streams) totals.push_back(s.totals);
        return totals;
    }

private:
    void Work() {
        std::unique_lock<std::mutex> lock(m_mutex);
        for (;;) {
            if (m_error || m_closed == m_streams.size()) break;
            Stream* stream = nullptr;
            if (m_in_flight < m_cap) {
                for (Stream* s : m_active) {
                    if (s->next_issue < s->round_end) {
                        stream = s;
                        break;
                    }
                }
                // all open samples wait (or, streamed, none is open): open the next
                if (!stream && !m_opening && m_next_sample < m_streams.size() && (!m_streaming || m_active.empty())) {
                    try {
                        Open(m_next_sample++, lock);
                    } catch (...) {
                        Fail(std::current_exception());
                    }
                    continue;
                }
            }
            if (!stream) {
                m_cv.wait(lock);
                continue;
            }
            std::size_t const k = stream->next_issue++;
            ++m_in_flight;
            Item const item = stream->items[k];
            std::shared_ptr<std::promise<std::shared_ptr<Contigs const>>> load;
            std::shared_future<std::shared_ptr<Contigs const>> loaded;
            fs::path fasta;
            if (item.genome >= 0) fasta = m_job.Fasta(stream->index, item.genome);
            if (!fasta.empty()) {
                auto [slot, fresh] = stream->genomes.try_emplace({item.round, item.genome});
                if (fresh) {
                    load = std::make_shared<std::promise<std::shared_ptr<Contigs const>>>();
                    slot->second.loaded = load->get_future().share();
                    slot->second.parts_left = item.parts;
                }
                loaded = slot->second.loaded;
            }
            lock.unlock();
            std::optional<Piece> piece;
            try {
                if (load) {
                    try {
                        load->set_value(LoadContigs(m_job.GenomeName(stream->index, item.genome), fasta,
                                                    m_job.MinContigLength()));
                    } catch (...) {
                        load->set_exception(std::current_exception());
                    }
                }
                piece = m_job.Make(stream->index, item, fasta.empty() ? nullptr : loaded.get().get());
            } catch (...) {
                lock.lock();
                Fail(std::current_exception());
                continue;
            }
            lock.lock();
            Complete(*stream, k, item, !fasta.empty(), std::move(*piece), lock);
        }
        m_cv.notify_all();
    }

    void Fail(std::exception_ptr error) {
        if (!m_error) m_error = error;
        m_cv.notify_all();
    }

    // Opens sample i's outputs (without the lock: a pipe waits for its reader) and plans its first round.
    void Open(std::size_t i, std::unique_lock<std::mutex>& lock) {
        Stream& s = m_streams[i];
        s.index = i;
        std::vector<Output> outputs;
        for (auto const& path : m_job.Outputs(i)) {
            Output out;
            out.path = path;
            out.packing = PackingOf(path);
            out.fifo = IsFifo(path);
            out.written = out.fifo ? path : fs::path(path.string() + ".partial");
            outputs.push_back(out);
        }
        m_opening = true;
        lock.unlock();
        std::exception_ptr error;
        for (auto& out : outputs) {
            try {
                if (!out.fifo && !out.path.parent_path().empty()) fs::create_directories(out.path.parent_path());
                out.file = std::fopen(out.written.c_str(), "wb");
                if (!out.file) throw std::runtime_error("cannot write " + out.written.string() + ": " + std::strerror(errno));
            } catch (...) {
                error = std::current_exception();
                break;
            }
        }
        lock.lock();
        m_opening = false;
        s.outputs = std::move(outputs);  // the destructor closes what was opened
        if (error) std::rethrow_exception(error);
        m_active.push_back(&s);
        PlanRound(s);
        Drain(s, lock);  // a sample without reads is done here
        m_cv.notify_all();
    }

    void PlanRound(Stream& s) {
        auto items = m_job.Plan(s.index, s.totals);
        if (items.empty()) {
            s.planned_all = true;
            return;
        }
        ++s.totals.rounds;
        for (auto& item : items) s.items.push_back(std::move(item));
        s.round_end = s.items.size();
        s.pieces.resize(s.items.size());
    }

    void Complete(Stream& s, std::size_t k, Item const& item, bool genome, Piece piece, std::unique_lock<std::mutex>& lock) {
        s.totals.reads += piece.reads;
        s.totals.template_bases += piece.template_bases;
        s.totals.read_bases += piece.read_bases;
        s.totals.errors += piece.errors;
        s.pieces[k] = std::move(piece);
        ++s.processed;
        if (genome) {
            auto slot = s.genomes.find({item.round, item.genome});
            if (slot != s.genomes.end() && --slot->second.parts_left == 0) s.genomes.erase(slot);
        }
        if (s.processed == s.round_end && !s.planned_all) {
            try {
                PlanRound(s);
            } catch (...) {
                Fail(std::current_exception());
            }
        }
        Drain(s, lock);
        m_cv.notify_all();
    }

    // Writes the sample's pieces that are next in order (each output's part of a piece, R1 before R2); one thread at
    // a time per sample, the lock released while it writes. Closes the files once the last piece is in.
    void Drain(Stream& s, std::unique_lock<std::mutex>& lock) {
        if (s.writing || s.closed) return;
        s.writing = true;
        while (!m_error && s.next_write < s.pieces.size() && s.pieces[s.next_write]) {
            Piece piece = std::move(*s.pieces[s.next_write]);
            s.pieces[s.next_write].reset();
            ++s.next_write;
            --m_in_flight;
            lock.unlock();
            std::string failed;
            for (std::size_t o = 0; o < s.outputs.size() && failed.empty(); ++o) {
                auto const& bytes = piece.bytes[o];
                if (std::fwrite(bytes.data(), 1, bytes.size(), s.outputs[o].file) != bytes.size()) {
                    failed = "writing " + s.outputs[o].written.string() + " failed: " + std::strerror(errno);
                }
            }
            lock.lock();
            if (!failed.empty()) Fail(std::make_exception_ptr(std::runtime_error(failed)));
        }
        if (!m_error && s.planned_all && s.next_write == s.items.size()) {
            lock.unlock();
            std::exception_ptr error;
            try {
                Close(s);
            } catch (...) {
                error = std::current_exception();
            }
            lock.lock();
            if (error) {
                Fail(error);
            } else {
                s.closed = true;
                ++m_closed;
                m_active.erase(std::find(m_active.begin(), m_active.end(), &s));
            }
        }
        s.writing = false;
    }

    void Close(Stream& s) {
        for (auto& out : s.outputs) {
            bool ok = true;
            if (out.packing == Packing::Bgzf) {
                ok = std::fwrite(protal::bgzf::kEof, 1, sizeof(protal::bgzf::kEof), out.file) == sizeof(protal::bgzf::kEof);
            } else if (s.items.empty()) {  // no reads: an empty frame, so that the file is valid zstd
                std::string const empty = Pack(std::string(), out.packing);
                ok = std::fwrite(empty.data(), 1, empty.size(), out.file) == empty.size();
            }
            ok = (std::fclose(out.file) == 0) && ok;
            out.file = nullptr;
            if (!ok) throw std::runtime_error("writing " + out.written.string() + " failed: " + std::strerror(errno));
            if (!out.fifo) fs::rename(out.written, out.path);
        }
    }

    Job& m_job;
    std::vector<Stream> m_streams;
    std::vector<Stream*> m_active;
    std::size_t m_next_sample = 0, m_closed = 0, m_in_flight = 0, m_cap = 4;
    int m_threads = 1;
    bool m_streaming = false, m_opening = false;
    std::exception_ptr m_error;
    std::mutex m_mutex;
    std::condition_variable m_cv;
};

}  // namespace

std::vector<Totals> Run(Job& job, int threads) {
    if (job.Samples() == 0) return {};
    Engine engine(job, threads);
    return engine.Run();
}

}  // namespace pipeline
}  // namespace protal::sim
