// ChunkRead with the overlap LongReadAligner asks for when the database's longest gene is L:
// L + 2 * (100 + L / 20). Counts genes of length L (every start, step 97) that do not lie wholly in
// the chunk whose core holds their centre, i.e. the "each gene in its owning chunk" rule.
#include <cstdio>
#include "LongReads.h"

using namespace protal;

int main() {
    for (int64_t L : { 4000, 20000, 30000, 31000, 35000, 50000 }) {
        size_t const requested = static_cast<size_t>(L + 2 * (100 + L / 20));
        for (size_t length : { 100000, 200000 }) {
            auto const chunks = ChunkRead(length, kMaxLongReadChunk, requested);
            size_t genes = 0, cut = 0;
            for (int64_t start = 0; start + L <= static_cast<int64_t>(length); start += 97) {
                ReadInterval const gene{ start, start + L };
                for (auto const& c : chunks) {
                    if (gene.TwiceCentre() >= 2 * c.core.start && gene.TwiceCentre() < 2 * c.core.end) {
                        genes++;
                        cut += gene.start < static_cast<int64_t>(c.offset) || gene.end > static_cast<int64_t>(c.offset + c.length);
                    }
                }
            }
            size_t const overlap = chunks.size() > 1 ? chunks[0].offset + chunks[0].length - chunks[1].offset : 0;
            std::printf("longest gene %6lld: overlap asked %6zu, read %6zu: %zu chunks, first overlap %6zu; genes %5zu, not wholly in their chunk %5zu\n",
                        static_cast<long long>(L), requested, length, chunks.size(), overlap, genes, cut);
        }
    }
}
