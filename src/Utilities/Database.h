// Database.h - where protal finds its database files: a directory of files (index.prx(.zst),
// reference.fna(.zst), reference.map, internal_taxonomy.dmp, unique_kmers.tsv, gene_conservation.tsv,
// suspect_copies.tsv and species_priors.tsv (optional), one presence model per read type: model_pe.xml, model_se.xml, model_PB.xml, model_ONT.xml,
// see ReadType.h), or the single-file database database.protal, which holds all of them as members.
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
//
// The directory frame holds its content uncompressed (RawFrame), so its size depends only on the
// members' names, not on their frame counts, and the database's writers put the models last: then
// --add_model replaces the models in place (ReplaceTail), writing the end of the file and the
// directory frame, instead of copying the ~20 GB before them (Write). Files written before this
// (a compressed directory, the models before gene_table.bin) are rewritten once by Write; readers
// of either kind read both.
#pragma once

#include "Zstd.h"

#include <algorithm>
#include <deque>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace protal::db {
    inline constexpr char kMagic[8] = {'P', 'R', 'O', 'T', 'A', 'L', 'D', 'B'};
    inline constexpr uint64_t kVersion = 1;
    inline const std::string kFileName = "database.protal";
    // Next to a database whose members ReplaceTail is replacing: the bytes it overwrites, until it is done.
    inline const std::string kJournalExtension = ".journal";

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

        // A zstd frame holding content as it is (raw blocks, no checksum), whose size depends on the content's
        // length alone (RawFrameSize).
        inline std::string RawFrame(std::string const& content) {
            std::string frame = {'\x28', '\xb5', '\x2f', '\xfd', '\xe0'};  // magic; 8-byte content size, single segment
            PutU64(frame, content.size());
            size_t constexpr kBlock = size_t{128} << 10;  // the largest block
            size_t pos = 0;
            do {
                size_t const n = std::min(kBlock, content.size() - pos);
                uint32_t const header = uint32_t(pos + n == content.size()) | uint32_t(n) << 3;  // last block?, raw, size
                for (int i = 0; i < 3; i++) frame.push_back(static_cast<char>((header >> (8 * i)) & 0xff));
                frame.append(content, pos, n);
                pos += n;
            } while (pos < content.size());
            return frame;
        }

        inline uint64_t RawFrameSize(uint64_t content) {
            uint64_t constexpr kBlock = uint64_t{128} << 10;
            return 13 + 3 * std::max<uint64_t>(1, (content + kBlock - 1) / kBlock) + content;
        }
    }

    inline bool SameFrames(zstd::SeekTable const& a, zstd::SeekTable const& b) {
        return std::equal(a.frames.begin(), a.frames.end(), b.frames.begin(), b.frames.end(), [](auto const& x, auto const& y) {
            return x.compressed_offset == y.compressed_offset && x.compressed_size == y.compressed_size &&
                   x.decompressed_offset == y.decompressed_offset && x.decompressed_size == y.decompressed_size;
        });
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
            bundle.m_directory = table->frames[0];
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
        // Frame 0, the directory.
        zstd::SeekTable::Frame const& DirectoryFrame() const { return m_directory; }

        Member const* Find(std::string const& name) const {
            for (auto const& member : m_members) {
                if (member.name == name) return &member;
            }
            return nullptr;
        }

    private:
        std::string m_path;
        zstd::SeekTable::Frame m_directory{};
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
        std::string missing;        // why nothing is at db ("does not exist", or why it cannot be accessed); empty if something is
    };

    inline Location Locate(std::string const& db) {
        namespace fs = std::filesystem;
        std::error_code ec;
        Location location;
        auto const status = fs::status(db, ec);
        if (fs::is_regular_file(status)) {
            auto const parent = fs::path(db).parent_path();
            location.dir = parent.empty() ? "." : parent.string();
            location.bundle = db;
            return location;
        }
        location.dir = db;
        if (status.type() == fs::file_type::not_found || !fs::status_known(status)) {  // unknown: e.g. no permission
            location.missing = fs::status_known(status) ? "does not exist" : "cannot be accessed: " + ec.message();
            return location;
        }
        std::string const bundle = (fs::path(db) / kFileName).string();
        bool const has_index = fs::exists(fs::path(db) / "index.prx", ec) || fs::exists(fs::path(db) / "index.prx.zst", ec);
        if (fs::exists(bundle, ec)) (has_index ? location.unused_bundle : location.bundle) = bundle;
        return location;
    }

    // A member whose frames are made in memory rather than read from a file: the index of protal --build, written
    // straight into the database (not into index.prx.zst first, then copied and both read back). Its number of frames
    // is known before (the directory frame goes first); write appends them to the database file; verify checks them in
    // it, where they are once written (`frames`: in the file, content offsets from the member's start); `done` follows
    // a successful check (--build frees the index there, before the other members are compressed). write and verify
    // return an error message, empty on success.
    struct Generated {
        uint64_t frames = 0;
        std::function<std::string(zstd::FrameWriter& out)> write;
        std::function<std::string(std::string const& path, zstd::SeekTable const& frames)> verify;
        std::function<void()> done;
    };

    // A member to write, from the file at path: a seekable zstd file's frames are copied as they are,
    // any other file (raw, or zstd without a seek table) is compressed into frames. With frames, only
    // those frames of path are copied (a member of another single-file database). With generated, the
    // member's frames come from it (path is only named in messages).
    struct Source {
        std::string name;
        std::string path;
        std::optional<zstd::SeekTable> frames = std::nullopt;
        std::shared_ptr<Generated> generated = nullptr;
    };

    namespace detail {
        struct Planned {
            Source source;
            std::optional<zstd::SeekTable> table;  // frames to copy, or none: compress
            uint64_t frames = 0;
            uint64_t size = 0;
            std::optional<zstd::SeekTable> written;  // a generated member's frames as written and checked
        };

        // What the database's members will be: per source its frames to copy (a seekable file's, or those given) or
        // the number it is compressed into, and its content size. nullopt with a message in error.
        inline std::optional<std::vector<Planned>> Plan(std::vector<Source> const& sources, zstd::Params const& params,
                                                        std::string& error) {
            if (params.frame_size == 0) {
                error = "a single-file database needs frames (frame size > 0)";
                return std::nullopt;
            }
            std::vector<Planned> plan;
            std::set<std::string> names;
            for (auto const& source : sources) {
                Planned p;
                p.source = source;
                if (!IsFileName(source.name) || !names.insert(source.name).second) {
                    error = "invalid or repeated member name '" + source.name + "'";
                    return std::nullopt;
                }
                if (source.generated) {  // its size once written
                    p.frames = source.generated->frames;
                    plan.push_back(std::move(p));
                    continue;
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
            return plan;
        }

        // The directory, frame 0's content, of the members `plan` lists.
        inline std::string Directory(std::vector<Planned> const& plan) {
            std::string directory(kMagic, kMagic + 8);
            PutU64(directory, kVersion);
            PutU64(directory, plan.size());
            uint64_t first = 1;
            for (auto const& p : plan) {
                PutU64(directory, p.source.name.size());
                directory += p.source.name;
                PutU64(directory, first);
                PutU64(directory, p.frames);
                first += p.frames;
            }
            return directory;
        }

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
                if (p.written) {  // generated: checked against its source when written, so its frames must be those
                    if (!SameFrames(member.frames, *p.written)) return what + "its frames are not those written";
                    continue;
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
        auto planned = detail::Plan(sources, params, error);
        if (!planned) return std::nullopt;
        auto& plan = *planned;
        std::string const directory = detail::Directory(plan);

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
            std::string const frame = detail::RawFrame(directory);
            if (!out.Add(frame.data(), frame.size(), directory.size())) return fail(out.Error());

            for (auto& p : plan) {
                uint64_t const before = out.Frames();
                if (p.source.generated) {
                    auto const& generated = *p.source.generated;
                    std::string const e = generated.write(out);
                    if (!e.empty()) return fail(p.source.name + ": " + e);
                    if (out.Frames() - before != p.frames) {
                        return fail(p.source.name + ": " + std::to_string(out.Frames() - before) + " frames written, " +
                                    std::to_string(p.frames) + " planned");
                    }
                    if (!out.Flush()) return fail(partial + ": " + out.Error());
                    auto const frames = out.Slice(before, out.Frames());
                    if (std::string const v = generated.verify(partial, frames); !v.empty()) {
                        return fail(p.source.name + " does not read back as written: " + v);
                    }
                    p.size = frames.DecompressedSize();
                    p.written = frames;
                    if (generated.done) generated.done();
                    continue;
                }
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
        std::filesystem::remove(path + kJournalExtension, ec);  // of the file just replaced, if any
        return written;
    }

    // Whether ReplaceTail can write sources into bundle in place, and from which member on: sources name the bundle's
    // members in the same order, those before that member as they are (the bundle's own frames), and the directory
    // frame keeps its size (raw, as Write writes it; a database written before is rewritten once). The members from
    // there on are held in memory, at most max_bytes of them. nullopt if not, or if nothing changes.
    inline std::optional<size_t> InPlaceFrom(Bundle const& bundle, std::vector<Source> const& sources,
                                             uint64_t max_bytes = uint64_t{1} << 30) {
        auto const& members = bundle.Members();
        if (sources.size() != members.size()) return std::nullopt;
        size_t first = 0;
        while (first < sources.size() && sources[first].name == members[first].name && sources[first].path == bundle.Path() &&
               sources[first].frames && SameFrames(*sources[first].frames, members[first].frames)) {
            first++;
        }
        if (first == sources.size()) return std::nullopt;
        uint64_t bytes = 0;
        for (size_t i = first; i < sources.size(); i++) {
            if (sources[i].name != members[i].name) return std::nullopt;
            if (sources[i].frames) {
                for (auto const& f : sources[i].frames->frames) bytes += f.compressed_size;
            } else {
                auto const size = zstd::UncompressedSize(sources[i].path);
                if (!size) return std::nullopt;
                bytes += *size;
            }
        }
        if (bytes > max_bytes) return std::nullopt;
        // The same names give a directory of the same length; raw, its frame then has the same size.
        auto const& directory = bundle.DirectoryFrame();
        if (directory.compressed_size != detail::RawFrameSize(directory.decompressed_size)) return std::nullopt;
        return first;
    }

    namespace detail {
        // What ReplaceTail overwrites, to put back: the directory frame and the file from tail_offset on; `check`, bytes
        // that ReplaceTail leaves as they are, before tail_offset, tells that a journal is of the file it is next to.
        struct Journal {
            static constexpr char kMagic[8] = {'P', 'R', 'O', 'T', 'A', 'L', 'J', '1'};
            uint64_t tail_offset = 0, check_offset = 0;
            std::string directory, check, tail;

            std::string Serialize() const {
                std::string out(kMagic, kMagic + 8);
                for (uint64_t v : {tail_offset, check_offset, uint64_t(directory.size()), uint64_t(check.size()), uint64_t(tail.size())}) {
                    PutU64(out, v);
                }
                return out + directory + check + tail;
            }

            static std::optional<Journal> Parse(std::vector<char> const& bytes) {
                if (bytes.size() < 48 || std::memcmp(bytes.data(), kMagic, 8) != 0) return std::nullopt;
                Journal j;
                size_t pos = 8;
                uint64_t sizes[3] = {};
                if (!GetU64(bytes, pos, j.tail_offset) || !GetU64(bytes, pos, j.check_offset)) return std::nullopt;
                for (auto& size : sizes) {
                    if (!GetU64(bytes, pos, size)) return std::nullopt;
                }
                if (sizes[0] + sizes[1] + sizes[2] != bytes.size() - pos || j.check_offset + sizes[1] != j.tail_offset ||
                    sizes[0] > j.check_offset) {
                    return std::nullopt;
                }
                std::pair<std::string*, uint64_t> const parts[] = {{&j.directory, sizes[0]}, {&j.check, sizes[1]}, {&j.tail, sizes[2]}};
                for (auto const& [part, size] : parts) {
                    part->assign(bytes.data() + pos, size);
                    pos += size;
                }
                return j;
            }
        };

        inline std::string ReadAt(int fd, uint64_t offset, uint64_t size, bool& ok) {
            std::string bytes(size, '\0');
            ok = ok && zstd::PreadAll(fd, bytes.data(), size, offset);
            return bytes;
        }

        inline bool WriteAt(int fd, std::string const& bytes, uint64_t offset) {
            char const* data = bytes.data();
            size_t size = bytes.size();
            while (size > 0) {
                ssize_t const n = ::pwrite(fd, data, size, static_cast<off_t>(offset));
                if (n < 0 && errno == EINTR) continue;
                if (n <= 0) return false;
                data += n;
                size -= static_cast<size_t>(n);
                offset += static_cast<uint64_t>(n);
            }
            return true;
        }

        // A file descriptor, closed with it.
        struct Descriptor {
            int fd;
            explicit Descriptor(int fd) : fd(fd) {}
            Descriptor(Descriptor const&) = delete;
            Descriptor& operator=(Descriptor const&) = delete;
            ~Descriptor() { if (fd >= 0) ::close(fd); }
        };

        // Writes `tail` at tail_offset into the file open as fd, ends the file there, and writes `directory` at its start,
        // durably (fsync reports what a network file system's close would). An error message, empty on success.
        inline std::string Splice(int fd, std::string const& directory, uint64_t tail_offset, std::string const& tail) {
            bool const ok = WriteAt(fd, tail, tail_offset) && ::ftruncate(fd, static_cast<off_t>(tail_offset + tail.size())) == 0 &&
                            WriteAt(fd, directory, 0) && ::fsync(fd) == 0;
            return ok ? "" : std::string("writing failed: ") + std::strerror(errno);
        }

        inline bool WriteFile(std::string const& path, std::string const& bytes) {
            int const fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC | O_CLOEXEC, 0644);
            if (fd < 0) return false;
            bool const ok = WriteAt(fd, bytes, 0) && ::fsync(fd) == 0;
            return ::close(fd) == 0 && ok;
        }
    }

    // Replaces bundle's members from `first` on (InPlaceFrom) with sources[first..] in place: their frames, compressed
    // params.threads at a time (or copied, read before anything is overwritten), and the new seek table overwrite the
    // file from the old frames of member `first` on, and the directory frame, of the same size, gets the new frame
    // counts. The members before `first` are neither read nor written. Until the result is checked (the new frames
    // read back as written, the compressed members' content as their sources), path.journal holds the old bytes; if
    // anything fails they are written back (the database is then unchanged), and a run stopped in between leaves
    // the journal, from which RestoreFromJournal writes them back. Returns the file's size, or nullopt with a message
    // in error.
    inline std::optional<uint64_t> ReplaceTail(Bundle const& bundle, std::vector<Source> const& sources, size_t first,
                                               zstd::Params const& params, std::string& error) {
        std::string const& path = bundle.Path();
        // Until the splice below, nothing in the file has changed.
        auto unchanged = [&error]() -> std::optional<uint64_t> {
            error += " (the database is unchanged)";
            return std::nullopt;
        };
        auto const planned = detail::Plan(sources, params, error);
        if (!planned) return unchanged();
        auto const& plan = *planned;
        auto const table = zstd::ReadSeekTable(path, error);
        if (!table) {
            if (error.empty()) error = path + " has no seek table";
            return unchanged();
        }
        uint64_t kept = 1;  // frames before member first's, with the directory
        for (size_t i = 0; i < first; i++) kept += plan[i].frames;
        auto const& last = table->frames.back();
        uint64_t const tail_offset = kept < table->frames.size() ? table->frames[kept].compressed_offset
                                                                 : last.compressed_offset + last.compressed_size;
        std::string const listing = detail::Directory(plan);
        std::string const directory = detail::RawFrame(listing);
        if (directory.size() != table->frames[0].compressed_size || listing.size() != table->frames[0].decompressed_size ||
            kept > table->frames.size()) {
            error = "the directory frame would change its size (the file needs rewriting)";
            return unchanged();
        }

        // The members' new frames: copied ones read now, the others' content read and compressed.
        struct Frame {
            std::string bytes;
            uint64_t content_size = 0;
        };
        std::vector<Frame> frames;
        std::deque<std::string> contents;
        std::vector<std::string const*> content_of(plan.size(), nullptr);
        std::vector<std::pair<size_t, std::string_view>> to_compress;  // frame, content
        for (size_t i = first; i < plan.size(); i++) {
            auto const& p = plan[i];
            size_t const before = frames.size();
            if (p.table) {
                int const fd = ::open(p.source.path.c_str(), O_RDONLY | O_CLOEXEC);
                if (fd < 0) {
                    error = "cannot open " + p.source.path + ": " + std::strerror(errno);
                    return unchanged();
                }
                bool ok = true;
                for (auto const& f : p.table->frames) frames.push_back({detail::ReadAt(fd, f.compressed_offset, f.compressed_size, ok), f.decompressed_size});
                ::close(fd);
                if (!ok) {
                    error = "read error in " + p.source.path;
                    return unchanged();
                }
            } else {
                zstd::InputFile in(p.source.path);
                std::string& content = contents.emplace_back();
                content_of[i] = &content;
                std::vector<char> buffer(size_t{1} << 20);
                try {
                    while (in.IsOpen() && in.Stream()) {
                        in.Stream().read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
                        content.append(buffer.data(), static_cast<size_t>(in.Stream().gcount()));
                    }
                } catch (std::exception const& e) {
                    error = "reading " + p.source.path + ": " + e.what();
                    return unchanged();
                }
                if (!in.IsOpen() || in.Stream().bad()) {
                    error = "cannot read " + p.source.path;
                    return unchanged();
                }
                for (uint64_t at = 0; at < content.size(); at += params.frame_size) {
                    auto const part = std::string_view(content).substr(at, params.frame_size);
                    to_compress.emplace_back(frames.size(), part);
                    frames.push_back({{}, part.size()});
                }
            }
            if (frames.size() - before != p.frames) {
                error = p.source.path + " changed while it was read";
                return unchanged();
            }
        }
        {
            using Context = std::unique_ptr<ZSTD_CCtx, size_t (*)(ZSTD_CCtx*)>;
            std::vector<Context> contexts;
            for (size_t w = 0; w < zstd::WorkerCount(to_compress.size(), params.threads); w++) {
                contexts.emplace_back(zstd::MakeCCtx(params, error), &ZSTD_freeCCtx);
                if (!contexts.back()) return unchanged();
            }
            error = zstd::ParallelFor(to_compress.size(), params.threads, [&](size_t u, size_t worker) -> std::string {
                auto const& [index, content] = to_compress[u];
                std::string& out = frames[index].bytes;
                out.resize(ZSTD_compressBound(content.size()));
                size_t const n = ZSTD_compress2(contexts[worker].get(), out.data(), out.size(), content.data(), content.size());
                if (ZSTD_isError(n)) return std::string("compressing: ") + ZSTD_getErrorName(n);
                out.resize(n);
                return "";
            });
            if (!error.empty()) return unchanged();
        }
        std::string entries, tail;
        for (uint64_t f = 0; f < kept; f++) {
            zstd::PutLE32(entries, static_cast<uint32_t>(table->frames[f].compressed_size));
            zstd::PutLE32(entries, static_cast<uint32_t>(table->frames[f].decompressed_size));
        }
        for (auto const& frame : frames) {
            if (frame.bytes.size() > 0xffffffffu || frame.content_size > 0xffffffffu) {
                error = "a frame is larger than 4 GB";
                return unchanged();
            }
            zstd::PutLE32(entries, static_cast<uint32_t>(frame.bytes.size()));
            zstd::PutLE32(entries, static_cast<uint32_t>(frame.content_size));
            tail += frame.bytes;
        }
        tail += zstd::SeekTableFrame(entries);

        // The old bytes into the journal, durably, before any is overwritten.
        detail::Descriptor const file(::open(path.c_str(), O_RDWR | O_CLOEXEC));
        if (file.fd < 0) {
            error = "cannot open " + path + " for writing: " + std::strerror(errno);
            return unchanged();
        }
        std::error_code ec;
        uint64_t const old_size = std::filesystem::file_size(path, ec);
        if (ec || old_size < tail_offset) {
            error = "cannot read the size of " + path;
            return unchanged();
        }
        detail::Journal journal;
        journal.tail_offset = tail_offset;
        journal.check_offset = std::max<uint64_t>(directory.size(), tail_offset >= 4096 ? tail_offset - 4096 : 0);
        bool read = true;
        journal.directory = detail::ReadAt(file.fd, 0, directory.size(), read);
        journal.check = detail::ReadAt(file.fd, journal.check_offset, tail_offset - journal.check_offset, read);
        journal.tail = detail::ReadAt(file.fd, tail_offset, old_size - tail_offset, read);
        if (!read) {
            error = "cannot read " + path;
            return unchanged();
        }
        std::string const journal_path = path + kJournalExtension;
        if (!detail::WriteFile(journal_path, journal.Serialize())) {
            std::filesystem::remove(journal_path, ec);
            error = "cannot write " + journal_path;
            return unchanged();
        }

        std::string problem = detail::Splice(file.fd, directory, tail_offset, tail);
        if (problem.empty()) {
            // The members before `first` are where they were; the others' frames read back as written, and the compressed
            // ones' content as their sources'.
            std::string open_error;
            auto const now = Bundle::Open(path, open_error);
            int const fd = file.fd;
            if (!now) problem = open_error.empty() ? "it is not a single-file database" : open_error;
            else if (now->Members().size() != plan.size()) problem = "it has " + std::to_string(now->Members().size()) + " members";
            size_t next = 0;
            for (size_t i = 0; problem.empty() && i < plan.size(); i++) {
                auto const& member = now->Members()[i];
                auto const& p = plan[i];
                std::string const what = "member " + member.name + ": ";
                if (member.name != p.source.name || member.frames.frames.size() != p.frames || member.Size() != p.size) {
                    problem = what + "does not match its source " + p.source.path;
                } else if (i < first) {
                    if (!SameFrames(member.frames, *p.table)) problem = what + "its frames moved";
                } else {
                    for (auto const& f : member.frames.frames) {
                        bool ok = true;
                        if (detail::ReadAt(fd, f.compressed_offset, f.compressed_size, ok) != frames[next++].bytes || !ok) {
                            problem = what + "a frame does not read back as written";
                            break;
                        }
                    }
                    if (problem.empty() && content_of[i]) {
                        std::string read_error;
                        auto const content = DbFile::InBundle(*now, member.name).ReadAll(read_error);
                        if (!content) problem = what + read_error;
                        else if (*content != *content_of[i]) problem = what + "its content differs from " + p.source.path;
                    }
                }
            }
        }
        if (!problem.empty()) {
            std::string const back = detail::Splice(file.fd, journal.directory, tail_offset, journal.tail);
            if (back.empty()) {
                std::filesystem::remove(journal_path, ec);
                error = problem + " (the old content is back: the database is unchanged)";
            } else {
                error = problem + "; writing the old content back failed too (" + back + "): " + journal_path + " holds it, and "
                        "protal --add_model with this database writes it back first";
            }
            return std::nullopt;
        }
        std::filesystem::remove(journal_path, ec);
        return tail_offset + tail.size();
    }

    // A database whose ReplaceTail stopped before it was done (its journal is next to it, and it may not open):
    // writes the journal's old bytes back, if the journal is of this file, and removes the journal. Returns false
    // with a message in error if it cannot.
    inline bool RestoreFromJournal(std::string const& path, std::string& error) {
        std::string const journal_path = path + kJournalExtension;
        std::vector<char> bytes;
        {
            std::ifstream is(journal_path, std::ios::binary);
            bytes.assign(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>());
        }
        auto const journal = detail::Journal::Parse(bytes);
        if (!journal) {
            error = journal_path + " is not a complete journal";
            return false;
        }
        {
            detail::Descriptor const file(::open(path.c_str(), O_RDWR | O_CLOEXEC));
            if (file.fd < 0) {
                error = "cannot open " + path + " for writing: " + std::strerror(errno);
                return false;
            }
            bool ok = true;
            if (detail::ReadAt(file.fd, journal->check_offset, journal->check.size(), ok) != journal->check || !ok) {
                error = journal_path + " is not of " + path + " (the bytes it keeps do not match)";
                return false;
            }
            error = detail::Splice(file.fd, journal->directory, journal->tail_offset, journal->tail);
            if (!error.empty()) return false;
        }
        std::string open_error;
        if (!Bundle::Open(path, open_error)) {
            error = "with the old bytes back, it still does not open" + (open_error.empty() ? "" : ": " + open_error);
            return false;
        }
        std::error_code ec;
        std::filesystem::remove(journal_path, ec);
        return true;
    }
}
