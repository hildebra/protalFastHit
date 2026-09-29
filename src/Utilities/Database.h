// Database.h - where protal finds its database files: a directory of files (index.prx(.zst),
// reference.fna(.zst), reference.map, internal_taxonomy.dmp, unique_kmers.tsv, one presence model per
// read type: model_pe.xml, model_se.xml, model_PB.xml, model_ONT.xml, see ReadType.h), or the
// single-file database database.protal, which holds all of them as members.
//
// database.protal is a seekable zstd file (Zstd.h):
//   frame 0       directory: "PROTALDB", version, member count, then per member its name, first
//                 frame and number of frames
//   frames 1..    each member's frames in turn: those of index.prx.zst (the index's column chunks,
//                 IndexCodec.h) and reference.fna.zst (64 MB frames) as these files hold them, the
//                 other files compressed the same way
//   seek table
// A member's frames are a slice of the seek table, so it is read like the .zst file it came from:
// the index and the reference load with -t threads as before. zstd -d on the whole file gives the
// members' contents one after another, which is of no use; protal --unpack_db writes the files.
#pragma once

#include "Zstd.h"

#include <filesystem>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace protal::db {
    inline constexpr char kMagic[8] = {'P', 'R', 'O', 'T', 'A', 'L', 'D', 'B'};
    inline constexpr uint64_t kVersion = 1;
    inline const std::string kFileName = "database.protal";

    struct Member {
        std::string name;
        zstd::SeekTable frames;  // offsets in the database file; content offsets from the member's start
        uint64_t Size() const { return frames.DecompressedSize(); }
    };

    // Member names are plain file names (--unpack_db writes each member as <dir>/<name>).
    inline bool IsFileName(std::string const& name) {
        return !name.empty() && name.size() <= 255 && name != "." && name != ".." &&
               name.find_first_of(std::string("/\\\0", 3)) == std::string::npos;
    }

    namespace detail {
        inline void PutU64(std::string& out, uint64_t v) {
            for (int i = 0; i < 8; i++) out.push_back(static_cast<char>((v >> (8 * i)) & 0xff));
        }

        inline bool GetU64(std::vector<char> const& data, size_t& pos, uint64_t& v) {
            if (pos + 8 > data.size()) return false;
            v = 0;
            for (int i = 0; i < 8; i++) v |= uint64_t(static_cast<unsigned char>(data[pos + i])) << (8 * i);
            pos += 8;
            return true;
        }

        // True if the (decompressed) file starts with the magic, e.g. a database.protal whose seek
        // table was cut off.
        inline bool StartsWithMagic(std::string const& path) {
            zstd::InputFile in(path);
            char magic[8] = {};
            if (!in.IsOpen()) return false;
            try {
                in.Stream().read(magic, 8);
            } catch (std::exception const&) {
                return false;
            }
            return in.Stream().gcount() == 8 && std::memcmp(magic, kMagic, 8) == 0;
        }
    }

    // An opened single-file database: its members and their frames.
    class Bundle {
    public:
        // The single-file database at path; nullopt with error empty if path is not one, with error
        // set if it is one but cannot be read.
        static std::optional<Bundle> Open(std::string const& path, std::string& error) {
            if (!zstd::IsCompressed(path)) return std::nullopt;
            auto const table = zstd::ReadSeekTable(path, error);
            if (!table) {
                if (error.empty() && detail::StartsWithMagic(path)) {
                    error = "the seek table at its end is missing (truncated or corrupt file?)";
                }
                return std::nullopt;
            }
            if (table->frames.empty()) return std::nullopt;
            // The directory is small; a large first frame is another file's (e.g. reference.fna.zst).
            if (table->frames[0].decompressed_size > (uint64_t{16} << 20)) {
                if (detail::StartsWithMagic(path)) error = "invalid directory: larger than 16 MB (corrupt file?)";
                return std::nullopt;
            }
            int const fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
            if (fd < 0) {
                error = std::string("cannot open the file: ") + std::strerror(errno);
                return std::nullopt;
            }
            ZSTD_DCtx* dctx = ZSTD_createDCtx();
            std::vector<char> input, data;
            std::string const read_error = dctx ? zstd::ReadFrame(fd, *table, 0, input, data, dctx)
                                                : "cannot allocate a zstd decompression context";
            if (dctx) ZSTD_freeDCtx(dctx);
            ::close(fd);
            if (!read_error.empty()) {
                // A frame 0 that is not the directory is another file (e.g. index.prx.zst) unless the
                // magic says otherwise.
                if (detail::StartsWithMagic(path)) error = "its directory cannot be read: " + read_error;
                return std::nullopt;
            }
            if (data.size() < 8 || std::memcmp(data.data(), kMagic, 8) != 0) return std::nullopt;

            auto fail = [&](std::string const& what) {
                error = "invalid directory: " + what + " (corrupt file?)";
                return std::nullopt;
            };
            size_t pos = 8;
            uint64_t version = 0, count = 0;
            if (!detail::GetU64(data, pos, version)) return fail("too short");
            if (version != kVersion) {
                error = "version " + std::to_string(version) + " of the single-file format; this protal reads version " +
                        std::to_string(kVersion) + " (a newer protal wrote it?)";
                return std::nullopt;
            }
            if (!detail::GetU64(data, pos, count) || count > table->frames.size()) return fail("member count");
            Bundle bundle;
            bundle.m_path = path;
            std::set<std::string> names;
            uint64_t next = 1;
            for (uint64_t i = 0; i < count; i++) {
                uint64_t length = 0, first = 0, frames = 0;
                if (!detail::GetU64(data, pos, length) || length == 0 || length > 4096 || pos + length > data.size()) {
                    return fail("member name");
                }
                Member member;
                member.name.assign(data.data() + pos, length);
                pos += length;
                if (!IsFileName(member.name) || member.name == kFileName) return fail("member name '" + member.name + "' is not a file name");
                if (!detail::GetU64(data, pos, first) || !detail::GetU64(data, pos, frames)) return fail("member frames");
                if (first != next || frames > table->frames.size() - first) return fail("the members do not tile the file");
                if (!names.insert(member.name).second) return fail("member " + member.name + " is listed twice");
                next += frames;
                uint64_t const base = frames > 0 ? table->frames[first].decompressed_offset : 0;
                for (uint64_t f = first; f < first + frames; f++) {
                    auto frame = table->frames[f];
                    frame.decompressed_offset -= base;
                    member.frames.frames.push_back(frame);
                }
                bundle.m_members.push_back(std::move(member));
            }
            if (next != table->frames.size()) return fail("the members do not cover the file");
            if (pos != data.size()) return fail("unexpected data after the member list");
            return bundle;
        }

        std::string const& Path() const { return m_path; }
        std::vector<Member> const& Members() const { return m_members; }

        Member const* Find(std::string const& name) const {
            for (auto const& member : m_members) {
                if (member.name == name) return &member;
            }
            return nullptr;
        }

    private:
        std::string m_path;
        std::vector<Member> m_members;
    };

    inline bool IsBundle(std::string const& path) {
        std::string error;
        return Bundle::Open(path, error).has_value();
    }

    // A database file: a file on disk (raw or .zst, as in a database directory), or a member of a
    // single-file database.
    class DbFile {
    public:
        DbFile() = default;

        // The file at path as it is (zstd::Resolve it first to find a .zst sibling).
        static DbFile OnDisk(std::string path) {
            DbFile file;
            file.m_path = path;
            file.m_name = std::move(path);
            file.m_exists = std::filesystem::exists(file.m_path);
            return file;
        }

        // Member `name` of bundle; Exists() is false if bundle has no such member.
        static DbFile InBundle(Bundle const& bundle, std::string const& name) {
            DbFile file;
            file.m_in_bundle = true;
            file.m_path = bundle.Path();
            file.m_name = name + " in " + bundle.Path();
            if (auto const* member = bundle.Find(name)) {
                file.m_exists = true;
                file.m_frames = member->frames;
            }
            return file;
        }

        bool Exists() const { return m_exists; }
        bool InBundle() const { return m_in_bundle; }
        // The file on disk: the file itself, or the single-file database.
        std::string const& Path() const { return m_path; }
        // For messages: the path, or "<member> in <database file>".
        std::string const& Name() const { return m_name; }
        // A member's frames in the database file.
        zstd::SeekTable const& Frames() const { return m_frames; }
        bool Compressed() const { return m_in_bundle || zstd::IsCompressed(m_path); }

        // Size of the content (uncompressed).
        std::optional<uint64_t> Size() const {
            if (!m_exists) return std::nullopt;
            return m_in_bundle ? std::optional<uint64_t>(m_frames.DecompressedSize()) : zstd::UncompressedSize(m_path);
        }

        // The content, read sequentially.
        std::unique_ptr<zstd::Input> Open() const {
            if (m_in_bundle) return std::make_unique<zstd::FramesInput>(m_path, m_frames, m_name);
            return std::make_unique<zstd::InputFile>(m_path);
        }

        // The content into sink with up to `threads` threads (zstd::ParallelRead).
        uint64_t ParallelRead(int threads, zstd::Sink& sink, std::string& error) const {
            if (m_in_bundle) return zstd::ParallelReadFrames(m_path, m_frames, threads, sink, error);
            return zstd::ParallelRead(m_path, threads, sink, error);
        }

        // The whole content (for small files); nullopt with a message in error.
        std::optional<std::string> ReadAll(std::string& error) const {
            if (!m_exists) {
                error = "does not exist";
                return std::nullopt;
            }
            auto input = Open();
            if (!input->IsOpen()) {
                error = "cannot open the file";
                return std::nullopt;
            }
            std::string content;
            std::vector<char> buffer(size_t{1} << 20);
            try {
                while (input->Stream()) {
                    input->Stream().read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                    content.append(buffer.data(), static_cast<size_t>(input->Stream().gcount()));
                }
            } catch (std::exception const& e) {
                error = e.what();
                return std::nullopt;
            }
            if (input->Stream().bad()) {
                error = "the file cannot be read or decompressed (truncated or corrupt file?)";
                return std::nullopt;
            }
            return content;
        }

    private:
        bool m_in_bundle = false;
        bool m_exists = false;
        std::string m_path, m_name;
        zstd::SeekTable m_frames;
    };

    // What --db names: a single-file database (a file, or a directory holding database.protal and no
    // index of its own), or a directory of database files.
    struct Location {
        std::string dir;            // the database directory; for a single file, the one it is in
        std::string bundle;         // the single-file database, empty for a directory of files
        std::string unused_bundle;  // a database.protal that the separate files next to it take precedence over
    };

    inline Location Locate(std::string const& db) {
        namespace fs = std::filesystem;
        std::error_code ec;
        Location location;
        if (fs::is_regular_file(db, ec)) {
            auto const parent = fs::path(db).parent_path();
            location.dir = parent.empty() ? "." : parent.string();
            location.bundle = db;
            return location;
        }
        location.dir = db;
        std::string const bundle = (fs::path(db) / kFileName).string();
        bool const has_index = fs::exists(fs::path(db) / "index.prx", ec) || fs::exists(fs::path(db) / "index.prx.zst", ec);
        if (fs::exists(bundle, ec)) (has_index ? location.unused_bundle : location.bundle) = bundle;
        return location;
    }

    // A member to write, from the file at path: a seekable zstd file's frames are copied as they are,
    // any other file (raw, or zstd without a seek table) is compressed into frames. With frames, only
    // those frames of path are copied (a member of another single-file database).
    struct Source {
        std::string name;
        std::string path;
        std::optional<zstd::SeekTable> frames = std::nullopt;
    };

    namespace detail {
        struct Planned {
            Source source;
            std::optional<zstd::SeekTable> table;  // frames to copy, or none: compress
            uint64_t frames = 0;
            uint64_t size = 0;
        };

        // Checks the database file at path against what was planned: the directory, then each
        // member: copied frames byte for byte, compressed members by content. Returns an error
        // message, empty if all matches.
        inline std::string Verify(std::string const& path, std::vector<Planned> const& plan, int threads) {
            std::string error;
            auto const bundle = Bundle::Open(path, error);
            if (!bundle) return error.empty() ? "it is not a single-file database" : error;
            if (bundle->Members().size() != plan.size()) return "it has " + std::to_string(bundle->Members().size()) + " members";
            std::vector<char> a(size_t{8} << 20), b(size_t{8} << 20);
            for (size_t i = 0; i < plan.size(); i++) {
                auto const& member = bundle->Members()[i];
                auto const& p = plan[i];
                std::string const what = "member " + member.name + ": ";
                if (member.name != p.source.name || member.frames.frames.size() != p.frames || member.Size() != p.size) {
                    return what + "does not match its source " + p.source.path;
                }
                if (p.table) {
                    int const src = ::open(p.source.path.c_str(), O_RDONLY | O_CLOEXEC);
                    int const dst = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
                    std::vector<std::vector<char>> ba(zstd::WorkerCount(p.frames, threads)), bb(ba.size());
                    std::string const e = src < 0 || dst < 0 ? std::string("cannot open the files")
                            : zstd::ParallelFor(p.frames, threads, [&](size_t f, size_t worker) -> std::string {
                        auto const& from = p.table->frames[f];
                        auto const& to = member.frames.frames[f];
                        if (from.compressed_size != to.compressed_size || from.decompressed_size != to.decompressed_size) {
                            return "frame " + std::to_string(f + 1) + " differs from " + p.source.path;
                        }
                        ba[worker].resize(from.compressed_size);
                        bb[worker].resize(to.compressed_size);
                        if (!zstd::PreadAll(src, ba[worker].data(), from.compressed_size, from.compressed_offset) ||
                            !zstd::PreadAll(dst, bb[worker].data(), to.compressed_size, to.compressed_offset)) {
                            return "read error in frame " + std::to_string(f + 1);
                        }
                        if (ba[worker] != bb[worker]) return "frame " + std::to_string(f + 1) + " differs from " + p.source.path;
                        return "";
                    });
                    if (src >= 0) ::close(src);
                    if (dst >= 0) ::close(dst);
                    if (!e.empty()) return what + e;
                    continue;
                }
                zstd::InputFile original(p.source.path);
                auto copy = DbFile::InBundle(*bundle, member.name).Open();
                if (!original.IsOpen() || !copy->IsOpen()) return what + "cannot open it or " + p.source.path;
                uint64_t compared = 0;
                try {
                    while (true) {
                        original.Stream().read(a.data(), static_cast<std::streamsize>(a.size()));
                        copy->Stream().read(b.data(), static_cast<std::streamsize>(b.size()));
                        std::streamsize const na = original.Stream().gcount(), nb = copy->Stream().gcount();
                        if (original.Stream().bad() || copy->Stream().bad()) return what + "a read failed";
                        if (na != nb || std::memcmp(a.data(), b.data(), static_cast<size_t>(na)) != 0) {
                            return what + "content differs from " + p.source.path + " after byte " + std::to_string(compared);
                        }
                        compared += static_cast<uint64_t>(na);
                        if (na == 0) break;
                    }
                } catch (std::exception const& e) {
                    return what + e.what();
                }
                if (compared != p.size) return what + std::to_string(compared) + " of " + std::to_string(p.size) + " bytes compared";
            }
            return "";
        }
    }

    // Writes a single-file database with the members `sources` (in this order) at path, via
    // path.partial, which is checked (detail::Verify) and then renamed. Members compressed here use
    // params (frame size, level, threads). Returns the bytes written, or nullopt with a message in
    // error.
    inline std::optional<uint64_t> Write(std::string const& path, std::vector<Source> const& sources,
                                         zstd::Params const& params, std::string& error) {
        if (params.frame_size == 0) {
            error = "a single-file database needs frames (frame size > 0)";
            return std::nullopt;
        }
        std::vector<detail::Planned> plan;
        std::set<std::string> names;
        for (auto const& source : sources) {
            detail::Planned p;
            p.source = source;
            if (!IsFileName(source.name) || !names.insert(source.name).second) {
                error = "invalid or repeated member name '" + source.name + "'";
                return std::nullopt;
            }
            if (!std::filesystem::exists(source.path)) {
                error = source.path + " does not exist";
                return std::nullopt;
            }
            if (source.frames) {
                p.table = source.frames;
            } else if (zstd::IsCompressed(source.path)) {
                p.table = zstd::ReadSeekTable(source.path, error);
                if (!error.empty()) {
                    error = source.path + ": " + error;
                    return std::nullopt;
                }
            }
            if (p.table) {
                p.frames = p.table->frames.size();
                p.size = p.table->DecompressedSize();
            } else {
                auto const size = zstd::UncompressedSize(source.path);
                if (!size) {
                    error = "cannot read " + source.path;
                    return std::nullopt;
                }
                p.size = *size;
                p.frames = (p.size + params.frame_size - 1) / params.frame_size;
            }
            plan.push_back(std::move(p));
        }

        std::string directory(kMagic, kMagic + 8);
        detail::PutU64(directory, kVersion);
        detail::PutU64(directory, plan.size());
        uint64_t first = 1;
        for (auto const& p : plan) {
            detail::PutU64(directory, p.source.name.size());
            directory += p.source.name;
            detail::PutU64(directory, first);
            detail::PutU64(directory, p.frames);
            first += p.frames;
        }

        std::string const partial = path + ".partial";
        std::error_code ec;
        auto fail = [&](std::string const& what) {
            error = what;
            std::filesystem::remove(partial, ec);
            return std::nullopt;
        };
        {
            zstd::FrameWriter out(partial);
            if (!out.Ok()) return fail(partial + ": " + out.Error());
            ZSTD_CCtx* cctx = zstd::MakeCCtx(params, error);
            if (!cctx) return fail(error);
            std::vector<char> frame(ZSTD_compressBound(directory.size()));
            size_t const n = ZSTD_compress2(cctx, frame.data(), frame.size(), directory.data(), directory.size());
            ZSTD_freeCCtx(cctx);
            if (ZSTD_isError(n)) return fail(std::string("compressing the directory: ") + ZSTD_getErrorName(n));
            if (!out.Add(frame.data(), n, directory.size())) return fail(out.Error());

            for (auto const& p : plan) {
                uint64_t const before = out.Frames();
                if (p.table) {
                    // Frames copied as they are, in large sequential reads.
                    int const fd = ::open(p.source.path.c_str(), O_RDONLY | O_CLOEXEC);
                    if (fd < 0) return fail("cannot open " + p.source.path + ": " + std::strerror(errno));
                    std::vector<char> buffer;
                    for (auto const& f : p.table->frames) {
                        buffer.resize(f.compressed_size);
                        if (!zstd::PreadAll(fd, buffer.data(), f.compressed_size, f.compressed_offset)) {
                            ::close(fd);
                            return fail("read error in " + p.source.path);
                        }
                        if (!out.Add(buffer.data(), buffer.size(), f.decompressed_size)) {
                            ::close(fd);
                            return fail(out.Error());
                        }
                    }
                    ::close(fd);
                } else {
                    zstd::InputFile in(p.source.path);
                    if (!in.IsOpen()) return fail("cannot open " + p.source.path);
                    zstd::StreamReader reader(in.Stream());
                    std::string compress_error;
                    try {
                        if (!zstd::CompressFramesTo(reader, out, params, compress_error)) {
                            return fail("compressing " + p.source.path + ": " + compress_error);
                        }
                    } catch (std::exception const& e) {
                        return fail("reading " + p.source.path + ": " + e.what());
                    }
                }
                if (out.Frames() - before != p.frames) {
                    return fail(p.source.path + " changed while it was written into " + partial);
                }
            }
            if (!out.Finish()) return fail(partial + ": " + out.Error());
        }
        std::string const verify_error = detail::Verify(partial, plan, params.threads);
        if (!verify_error.empty()) return fail("checking " + partial + ": " + verify_error);
        uint64_t const written = std::filesystem::file_size(partial, ec);
        std::filesystem::rename(partial, path, ec);
        if (ec) return fail("cannot rename " + partial + " to " + path + ": " + ec.message());
        return written;
    }
}
