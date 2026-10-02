// copybench - experiment only: how fast can the records file be appended behind a header? Copies SRC to the
// end of a new DEST that starts with HEADER_BYTES bytes, as SamOutput::AppendRecords does, on T threads that
// each take a range of SRC: "cfr" with copy_file_range at explicit offsets, "rw" with pread/pwrite through a
// 4 MB buffer. Prints the seconds and GB/s; DEST is removed afterwards unless KEEP is set.
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>

int main(int argc, char** argv) {
    if (argc < 6) { std::fprintf(stderr, "usage: copybench SRC DEST cfr|rw THREADS HEADER_BYTES\n"); return 2; }
    std::string const mode = argv[3];
    int const threads = std::atoi(argv[4]);
    off_t const head = std::atoll(argv[5]);
    int const in = ::open(argv[1], O_RDONLY);
    struct stat st{};
    if (in < 0 || ::fstat(in, &st) != 0) { std::perror("src"); return 1; }
    off_t const size = st.st_size;
    auto const t0 = std::chrono::steady_clock::now();
    int const out = ::open(argv[2], O_WRONLY | O_CREAT | O_TRUNC, 0666);
    if (out < 0) { std::perror("dest"); return 1; }
    std::vector<char> header(static_cast<size_t>(head), 'H');
    if (head > 0 && ::write(out, header.data(), header.size()) != static_cast<ssize_t>(header.size())) { std::perror("header"); return 1; }
    std::atomic<bool> failed{ false };
    off_t const step = (size + threads - 1) / threads;
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; t++) {
        pool.emplace_back([&, t] {
            off_t from = std::min<off_t>(size, step * t), to = std::min<off_t>(size, from + step);
            if (mode == "cfr") {
                loff_t in_off = from, out_off = head + from;
                while (in_off < to) {
                    ssize_t const n = ::copy_file_range(in, &in_off, out, &out_off, static_cast<size_t>(to - in_off), 0);
                    if (n <= 0) { if (n < 0 && errno == EINTR) continue; failed = true; return; }
                }
            } else {
                std::vector<char> buffer(size_t{4} << 20);
                for (off_t at = from; at < to;) {
                    ssize_t const n = ::pread(in, buffer.data(), std::min<off_t>(buffer.size(), to - at), at);
                    if (n <= 0) { failed = true; return; }
                    for (ssize_t done = 0; done < n;) {
                        ssize_t const w = ::pwrite(out, buffer.data() + done, n - done, head + at + done);
                        if (w <= 0) { failed = true; return; }
                        done += w;
                    }
                    at += n;
                }
            }
        });
    }
    for (auto& thread : pool) thread.join();
    ::close(out);
    double const s = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("mode\t%s\tthreads\t%d\tbytes\t%lld\tseconds\t%.3f\tGB_per_s\t%.2f\tfailed\t%d\n", mode.c_str(), threads,
                static_cast<long long>(size), s, size / s / 1e9, static_cast<int>(failed));
    if (!std::getenv("KEEP")) std::remove(argv[2]);
    return failed ? 1 : 0;
}
