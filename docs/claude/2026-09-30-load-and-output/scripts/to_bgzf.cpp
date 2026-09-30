// to_bgzf SRC DST THREADS: SRC (plain) written as BGZF with protal's bgzf::CompressFile.
#include <cstdio>
#include <cstdlib>
#include "Bgzf.h"

int main(int argc, char** argv) {
    if (argc != 4) { std::fprintf(stderr, "usage: to_bgzf SRC DST THREADS\n"); return 2; }
    std::string const error = protal::bgzf::CompressFile(argv[1], argv[2], std::atoi(argv[3]));
    if (!error.empty()) { std::fprintf(stderr, "%s\n", error.c_str()); return 1; }
    return 0;
}
