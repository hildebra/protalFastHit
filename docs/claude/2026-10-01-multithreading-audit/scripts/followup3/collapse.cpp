// collapse - experiment only: what FALLOC_FL_COLLAPSE_RANGE costs on a freshly written file. Writes PATH as
// 64 KB of room (a hole) and DATA_BYTES of data behind it, then collapses 60 KB of the room; prints the seconds of
// the write and of the collapse. With "sync" first, the data is written to disk before the collapse.
#include <fcntl.h>
#include <unistd.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: collapse PATH DATA_BYTES [sync]\n"); return 2; }
    size_t const bytes = std::strtoull(argv[2], nullptr, 10);
    bool const sync_first = argc > 3 && std::string(argv[3]) == "sync";
    auto const t0 = std::chrono::steady_clock::now();
    int const fd = ::open(argv[1], O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (fd < 0 || ::lseek(fd, 65536, SEEK_SET) < 0) { std::perror("open"); return 1; }
    std::vector<char> buffer(size_t{1} << 20, 'x');
    for (size_t done = 0; done < bytes;) {
        size_t const n = std::min(buffer.size(), bytes - done);
        if (::write(fd, buffer.data(), n) != static_cast<ssize_t>(n)) { std::perror("write"); return 1; }
        done += n;
    }
    if (sync_first) ::fdatasync(fd);
    auto const t1 = std::chrono::steady_clock::now();
    int const r = ::fallocate(fd, FALLOC_FL_COLLAPSE_RANGE, 4096, 65536 - 4096);
    auto const t2 = std::chrono::steady_clock::now();
    ::close(fd);
    std::printf("bytes\t%zu\tsync_first\t%d\twrite_s\t%.3f\tcollapse_s\t%.3f\tresult\t%s\n", bytes, sync_first,
                std::chrono::duration<double>(t1 - t0).count(), std::chrono::duration<double>(t2 - t1).count(),
                r == 0 ? "ok" : std::strerror(errno));
    std::remove(argv[1]);
    return 0;
}
