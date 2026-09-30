// ASan fuzz of index_codec::detail::DecodeChunk on mutated/truncated encoded chunks.
// Destination buffers are sized EXACTLY to the chunk (c.blocks*cpb keymap cells, c.values values),
// so any out-of-bounds read/write on corrupt input is caught by AddressSanitizer.
#include "Hash/IndexCodec.h"
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>
using namespace protal;
using namespace protal::index_codec;

struct SmallIndex {
    Layout layout;
    std::vector<uint16_t> keymap;
    std::vector<uint64_t> values;
    SmallIndex(uint64_t blocks, unsigned seed, double empty_share) {
        std::mt19937_64 rng(seed);
        layout.blocks = blocks;
        keymap.assign(layout.KeymapCells(), 0);
        uint64_t v = 0;
        for (uint64_t b = 0; b < blocks; b++) {
            uint16_t* cells = keymap.data() + b * layout.CellsPerBlock();
            std::memcpy(cells, &v, 8);
            uint64_t local = 0;
            for (uint64_t j = 0; j < layout.keys_per_block; j++) {
                cells[4 + j] = static_cast<uint16_t>(local);
                std::uniform_real_distribution<double> u(0, 1);
                uint64_t const c = u(rng) < empty_share ? 0 : 1 + rng() % (rng() % 4 == 0 ? 40 : 3);
                uint64_t const n = c >= 2 ? c + (c + 1) / 2 : c;
                uint64_t const taxid = 1 + rng() % 5000;
                for (uint64_t i = 0; i < n - c; i++) values.push_back(rng());
                for (uint64_t i = 0; i < c; i++) {
                    uint64_t const t = taxid + rng() % 20, g = 1 + rng() % 168, p = rng() % 4000, f = rng() % 4;
                    values.push_back(t << 40 | g << 20 | p | f << 60);
                }
                local += n;
            }
            v += local;
        }
        std::memcpy(keymap.data() + blocks * layout.CellsPerBlock(), &v, 8);
        layout.values = values.size();
    }
};

static long long ok_count = 0, err_count = 0;

static void TryDecode(std::vector<char> const& buf, Layout const& l, Chunk const& c) {
    uint64_t const cpb = l.CellsPerBlock();
    std::vector<uint16_t> km(c.blocks * cpb, 0xabcd);
    std::vector<uint64_t> vals(c.values, 0xabcdefULL);
    std::string e = detail::DecodeChunk(buf.data(), buf.size(), l, c, km.data(), vals.data());
    if (e.empty()) ok_count++; else err_count++;
}

int main() {
    std::mt19937_64 rng(12345);
    for (unsigned seed : {1u, 2u, 3u, 7u}) {
        for (double share : {0.2, 0.5, 0.9}) {
            SmallIndex idx(400, seed, share);
            std::string cerr_;
            auto chunks = detail::MakeChunks(idx.layout, idx.keymap.data(), 4096, cerr_);
            if (!chunks) { std::printf("MakeChunks failed: %s\n", cerr_.c_str()); return 1; }
            for (auto const& c : *chunks) {
                std::vector<char> good;
                detail::EncodeChunk(idx.layout, c, idx.keymap.data(), idx.values.data(), good);
                // sanity: the good encoding decodes
                TryDecode(good, idx.layout, c);
                // 1) single-byte flips across the header + a prefix of the body
                size_t const span = std::min<size_t>(good.size(), 3000);
                for (size_t i = 0; i < span; i++) {
                    for (unsigned char mask : {0x01, 0x80, 0xFF}) {
                        std::vector<char> m = good; m[i] ^= mask; TryDecode(m, idx.layout, c);
                    }
                }
                // 2) targeted header-field corruption (bytes 0..31): widths, entries/flex/filled
                for (int t = 0; t < 400; t++) {
                    std::vector<char> m = good;
                    m[0] = (char)(rng() % 4);            // mode 0..3
                    m[1] = (char)(rng() % 256);          // tw
                    m[2] = (char)(rng() % 256);          // gw
                    m[3] = (char)(rng() % 256);          // pw
                    uint64_t e = rng(), f = rng(), fi = rng();
                    std::memcpy(m.data() + 8, &e, 8);
                    std::memcpy(m.data() + 16, &f, 8);
                    std::memcpy(m.data() + 24, &fi, 8);
                    TryDecode(m, idx.layout, c);
                }
                // 3) truncations at every length
                for (size_t len = 0; len <= good.size(); len += (good.size() > 200 ? 7 : 1)) {
                    std::vector<char> m(good.begin(), good.begin() + len);
                    TryDecode(m, idx.layout, c);
                }
                // 4) random multi-byte smashes
                for (int t = 0; t < 300; t++) {
                    std::vector<char> m = good;
                    int nmut = 1 + rng() % 8;
                    for (int k = 0; k < nmut; k++) m[rng() % m.size()] ^= (char)(1 + rng() % 255);
                    TryDecode(m, idx.layout, c);
                }
            }
        }
    }
    std::printf("done: decoded_ok=%lld returned_error=%lld (no ASan abort => no OOB)\n", ok_count, err_count);
    return 0;
}
