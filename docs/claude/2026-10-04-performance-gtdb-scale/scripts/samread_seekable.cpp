// Compresses a plain SAM into protal's seekable .sam.zst layout (1 MB frames, level 3).
#include "Utilities/Zstd.h"
#include <iostream>
int main(int argc, char** argv) {
    if (argc != 3) { std::cerr << "usage: seekable in.sam out.sam.zst\n"; return 1; }
    std::string error;
    protal::zstd::Params params{ 3, 0, 6, uint64_t{1} << 20 };
    if (!protal::zstd::CompressFile(argv[1], argv[2], params, false, error)) { std::cerr << error << "\n"; return 1; }
    return 0;
}
